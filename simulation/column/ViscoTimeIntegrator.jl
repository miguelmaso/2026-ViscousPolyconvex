using Revise
using HyperFEM, HyperFEM.ComputationalModels.PostMetrics
using Gridap, Gridap.FESpaces, Gridap.Geometry
using GridapSolvers, GridapSolvers.NonlinearSolvers
using Printf
using Plots
using JLD2

default(palette = :seaborn_colorblind)


## Domain

function generate_tessellation(; width, height, hrefinement, xdivisions, args...)
  domain = (-0.5width, 0.5width, 0.0, height, -0.5width, 0.5width)
  ydivisions = round(Int, xdivisions * height / width)
  partition = hrefinement .* (xdivisions, ydivisions, xdivisions)
  geometry = CartesianDiscreteModel(domain, partition)
  labels = get_face_labeling(geometry)
  add_tag_from_tags!(labels, "bottom", CartesianTags.faceX0Z⁺)
  geometry
end


## Constitutive model

function build_model(; ν, args...)

  # VHB4910, calibrated parameters, see calibration/Calibration.jl (Random.seed!(1234))

  # Equilibrium branch
  μe = 13446.1  # [Pa]

  # Non-equilibrium branches
  μ1 = 11053.7  # [Pa]
  τ1 = 151.68   # [s]
  μ2 = 30480.1  # [Pa]
  τ2 = 7.0731   # [s]

  # Bulk modulus from the instantaneous shear modulus and the Poisson ratio
  μ0 = μe + μ1 + μ2
  κr = 2μ0 * (1 + ν) / (3 * (1 - 2ν))  # [Pa]

  equilibrium = NeoHookean3D(μ=μe, λ=0.0) + VolumetricEnergy(λ=κr)
  branch_1 = ViscousPolyconvex(μ=μ1, τ=τ1)
  branch_2 = ViscousPolyconvex(μ=μ2, τ=τ2)
  visco_model = GeneralizedMaxwell(equilibrium, branch_1, branch_2)
  return visco_model
end


## FEM solver

function solve_problem(data)

  pname = stem(@__FILE__)
  folder = abspath(dirname(@__FILE__), "results")
  outpath = joinpath(folder, pname)
  setupfolder(folder; remove=".vtu")

  model = build_model(; data...)

  F, H, J = get_Kinematics(Kinematics(Mechano, Solid))
  E       = get_Kinematics(Kinematics(Electro, Solid))

  geometry = generate_tessellation(; data...)

  writevtk(geometry, outpath * "_geometry")

  # Discrete domain, test FE

  Δt = data.Δt
  t_end = data.t_end
  order = data.prefinement
  degree = 2 * order
  Ω = Triangulation(geometry)
  dΩ = Measure(Ω, degree)

  solver = FESolver(NewtonSolver(LUSolver(); maxiter=20, atol=1e-8,  rtol=1e-8,  verbose=true))
  reffe = ReferenceFE(lagrangian, VectorValue{3,Float64}, order)
  Vu = TestFESpace(geometry, reffe, data.dirichlet_u, conformity=:H1)

  println("======================================\n")
  println("Degrees of freedom : $(Vu.nfree)\n")
  println("======================================")

  # Trial FE spaces and state variables

  Uu  = TrialFESpace(Vu, data.dirichlet_u)
  Uu⁻ = TrialFESpace(Vu, data.dirichlet_u)
  Uυ  = TrialFESpace(Vu, data.dirichlet_u)
  uh⁺ = FEFunction(Uu,  zero_free_values(Uu))
  uh⁻ = FEFunction(Uu⁻, zero_free_values(Uu))
  υh  = interpolate(data.velocity₀, Uυ)

  Fh  = F∘∇(uh⁺)'
  Fh⁻ = F∘∇(uh⁻)'
  Ah  = CellState(model, dΩ)
  Ph⁻ = CellState(zero(TensorValue{3,3,Float64}), dΩ)

  # Weak forms: residual and jacobian

  update_time_step!(model, Δt)
  Ψ, ∂Ψ∂F, ∂∂Ψ∂FF = model()
  D = Dissipation(model)
  P_updater(P, F, Fn, A...) = (true, ∂Ψ∂F(F, Fn, A...))

  function residual(u, v)
    ρ₀ = data.ρ₀
    a_mid = 2/Δt^2 * (u - uh⁻) - 2/Δt*υh
    P_mid = 1/2*(∂Ψ∂F ∘ (F ∘ ∇(u)', Fh⁻, Ah...) + Ph⁻)
    ∫( ρ₀ * a_mid · v + ∇(v)' ⊙ P_mid )dΩ
  end

  function jacobian(u, du, v)
    ρ₀ = data.ρ₀
    Da = 2/Δt^2 * du
    DP = 1/2*(∂∂Ψ∂FF ∘ (F ∘ ∇(u)', Fh⁻, Ah...)) ⊙ ∇(du)'
    ∫( ρ₀ * Da · v + ∇(v)' ⊙ DP )dΩ
  end

  # Post-processor

  fields = (:time, :KE, :EE, :Dvis, :Wext)
  metrics = NamedTuple(f => Float64[] for f in fields)

  function post_metrics!(m, step, time)
    ρ₀ = data.ρ₀
    Dvis⁻ = step > 0 ? m.Dvis[end] : 0.0
    Wext = sum(∫( 0.5ρ₀*υh·υh )dΩ)  # The column is clamped: the only energy input is the initial kinetic energy
    Wext⁻ = step > 0 ? m.Wext[end] : Wext
    υh⁺ = 2/Δt*(uh⁺ - uh⁻) - υh
    push!(m.time, time)
    push!(m.KE,   sum(∫( 0.5ρ₀*υh⁺·υh⁺ )dΩ))
    push!(m.EE,   sum(∫( Ψ∘(Fh, Fh⁻, Ah...) )dΩ))
    push!(m.Dvis, sum(∫( D∘(Fh, Fh⁻, Ah...) )dΩ) * Δt + Dvis⁻)
    push!(m.Wext, Wext⁻)
  end

  function post_vtk!(pvd, step, time)
    if mod(step, 10) == 0
      Ph = interpolate_L2_field(∂Ψ∂F ∘ (Fh, Fh⁻, Ah...), Ω, dΩ)
      Jh = interpolate_L2_field(J∘Fh, Ω, dΩ)
      pvd[time] = createvtk(Ω, outpath * @sprintf("_%03d", step), cellfields=["u" => uh⁺, "v" => υh, "J" => Jh, "P" => Ph], nsubcells=data.prefinement)
    end
  end

  # Time integration

  update_state!(P_updater, Ph⁻, Fh, Fh⁻, Ah...)
  update_state!(model, Ah, Fh, Fh⁻)

  createpvd(outpath) do pvd
    step = 0
    time = 0.0
    post_vtk!(pvd, step, time)
    post_metrics!(metrics, step, time)
    try
      while time < t_end - 1e-3Δt
        step += 1
        time += Δt
        printstyled(@sprintf("Step: %i\nTime: %.4f s\n", step, time), color=:green, bold=true)

        #-----------------------------------------
        # Solve (homogeneous Dirichlet at the bottom)
        #-----------------------------------------
        op = FEOperator(residual, jacobian, Uu, Vu)
        solve!(uh⁺, solver, op)

        #-----------------------------------------
        # Post processing
        #-----------------------------------------
        post_vtk!(pvd, step, time)
        post_metrics!(metrics, step, time)

        #-----------------------------------------
        # Update old step
        #-----------------------------------------
        update_state!(P_updater, Ph⁻, Fh, Fh⁻, Ah...)
        update_state!(model, Ah, Fh, Fh⁻)
        update_velocity!(υh, uh⁺, uh⁻, Δt)
        update_displacements!(uh⁻, uh⁺)
      end
    catch e
      rethrow(e)
    finally
      @save "$(outpath)_metrics.jld2" metrics  # Save the time evolution
    end
  end
  return (; metrics, uh⁺)
end


## Problem definition & solve
# Twisting column: 1 x 6 x 1 cm^3, clamped at Y = 0, free at the top,
# initial velocity v₀ = 1100 sin(πY/0.06) [Z, 0, -X] m/s

problem_data = let
  width = 0.01   # [m]
  height = 0.06  # [m]
  ν = 0.49
  ρ₀ = 960.0     # [kg/m^3]
  xdivisions = 2
  hrefinement = 1
  prefinement = 2
  t_end = 0.2    # [s]
  Δt = 1e-3      # [s]

  ω₀ = 1100.0    # [1/s]
  velocity₀(X) = ω₀ * sin(π * X[2] / height) * VectorValue(X[3], 0.0, -X[1])

  dirichlet_u = DirichletBC(["bottom"], [[0.0, 0.0, 0.0]])

  (; width, height, ν, ρ₀, hrefinement, xdivisions, prefinement, t_end, Δt, velocity₀, dirichlet_u)
end

metrics, uh = solve_problem(problem_data)

Wdiff = metrics.Dvis + metrics.KE + metrics.EE - metrics.Wext

@show Wdiff[1]
@show sum(Wdiff) * problem_data.Δt

let p = areaplot(metrics.time, [metrics.Dvis metrics.KE metrics.EE], label=["Dissipation" "Kinetic" "Elastic"], alpha=0.9)
  plot!(metrics.time, metrics.Wext, label="External", color=:black, lw=3, style=:dot)
  plot!(xlabel="Time [s]", ylabel="Energy [J]", legend=:outertop, legend_column=-1)
  display(p)
end
