using Revise
using HyperFEM, HyperFEM.ComputationalModels.PostMetrics
using Gridap, Gridap.FESpaces, Gridap.Geometry
using GridapSolvers, GridapSolvers.NonlinearSolvers
using Printf
using Plots
using JLD2

default(palette = :seaborn_colorblind)


## Discrete domain

function generate_tessellation(; width, long, hrefinement, args...)
  domain = (-0.5width, 0.5width, 0.0, long, -0.5width, 0.5width)
  partition = hrefinement .* (1, round(long/width), 1)
  geometry = CartesianDiscreteModel(domain, partition)
  labels = get_face_labeling(geometry)
  add_tag_from_tags!(labels, "top",    CartesianTags.faceX1Z⁺)
  add_tag_from_tags!(labels, "bottom", CartesianTags.faceX0Z⁺)
  geometry
end


## Constitutive model

function build_model(; args...)

  # Equilibrium branch
  μe = 13.4e3  # [Pa]

  # Non-equilibrium branches
  μ1 = 31.2e3  # [Pa]
  τ1 = 6.45    # [s]
  μ2 = 11.6e3  # [Pa]
  τ2 = 144.1   # [s]

  # Bulk modulus
  bulk_modulus(μ, ν) = 2/3*(1 - ν)/(1 - 2ν)*μ
  κr = bulk_modulus(μe + μ1 + μ2, 0.49)
  # κr = 91.43e3 # [Pa]  ## 2.5e6   # [Pa]

  equilibrium = NeoHookean3D(μ=μe, λ=0.0) + VolumetricEnergy(λ=κr)
  branch_1 = ViscousPolyconvex(μ=μ1, τ=τ1)
  branch_2 = ViscousPolyconvex(μ=μ2, τ=τ2)
  viscous_model = GeneralizedMaxwell(equilibrium, branch_1, branch_2)
  return viscous_model
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
  Vu = TestFESpace(geometry, reffe, data.dirichlet, conformity=:H1)

  println("======================================\n")
  println("Degrees of freedom : $(Vu.nfree)\n")
  println("======================================")

  # Trial FE spaces and state variables

  Uu  = TrialFESpace(Vu, data.dirichlet)
  Uu⁻ = TrialFESpace(Vu, data.dirichlet)
  Uυ  = TrialFESpace(Vu, get_time_derivative(data.dirichlet, 0.0))
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

  function residual(u, v)  # TODO: convertir en una función del tiempo/Λ
    ρ₀ = 960.0
    a_mid = 2/Δt^2 * (u - uh⁻) - 2/Δt*υh
    P_mid = 1/2*(∂Ψ∂F ∘ (F ∘ ∇(u)', Fh⁻, Ah...) + Ph⁻)
    ∫( ρ₀ * a_mid · v + ∇(v)' ⊙ P_mid )dΩ
  end

  function jacobian(u, du, v)
    ρ₀ = 960.0
    Da = 2/Δt^2 * du
    DP = 1/2*(∂∂Ψ∂FF ∘ (F ∘ ∇(u)', Fh⁻, Ah...)) ⊙ ∇(du)'
    ∫( ρ₀ * Da · v + ∇(v)' ⊙ DP )dΩ
  end

  # Post-processor

  fields = (:time, :KE, :EE, :Dvis, :Wext)
  metrics = NamedTuple(f => Float64[] for f in fields)

  function post_metrics!(data, step, time)
    ρ₀ = 960.0
    b = assemble_vector(v -> residual(uh⁺, v), DirichletFESpace(Vu))[:]
    Δu_dir = get_dirichlet_dof_values(Uu) - get_dirichlet_dof_values(Uu⁻)
    Dvis⁻ = step > 0 ? data.Dvis[end] : 0.0
    Wext⁻ = step > 0 ? data.Wext[end] : sum(∫( 0.5ρ₀*υh·υh )dΩ)
    υh⁺ = 2/Δt*(uh⁺ - uh⁻) - υh
    push!(data.time, time)
    push!(data.KE,   sum(∫( 0.5ρ₀*υh⁺·υh⁺ )dΩ))
    push!(data.EE,   sum(∫( Ψ∘(Fh, Fh⁻, Ah...) )dΩ))
    push!(data.Dvis, sum(∫( D∘(Fh, Fh⁻, Ah...) )dΩ) * Δt + Dvis⁻)
    push!(data.Wext, b · Δu_dir + Wext⁻)
  end

  function post_vtk!(pvd, step, time)
    if mod(step, 5) == 0
      P = ∂Ψ∂F ∘ (Fh, Fh⁻, Ah...)
      Ph = interpolate_L2_field(P, dΩ)
      Jh = interpolate_L2_field(J∘Fh, dΩ)
      ph = interpolate_L2_field(1/(3*(J∘Fh)) * P ⊙ Fh, dΩ)
      pvd[time] = createvtk(Ω, outpath * @sprintf("_%03d", step), cellfields=["u" => uh⁺, "v" => υh, "J" => Jh, "P" => Ph, "p" => ph], nsubcells=data.prefinement)
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
      while time < t_end
        step += 1
        time += Δt
        printstyled(@sprintf("Step: %i\nTime: %.4f s\n", step, time), color=:green, bold=true)

        #-----------------------------------------
        # Apply boundary conditions and solve
        #-----------------------------------------
        TrialFESpace!(Uu, data.dirichlet, time)

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

problem_data = let
  width = 0.01  # [m]
  long = 0.06   # [m]
  ω₀ = 1100.0   # [1/s]
  velocity₀((X,Y,Z),) = ω₀ * sin(π * Y / long) * VectorValue(Z, 0.0, -X)
  hrefinement = 2
  prefinement = 2
  t_end = 0.1
  linear_velocity = 0.5width * √2 * ω₀
  Δt = 0.2 * width / (prefinement * hrefinement * linear_velocity)
  @show Δt

  dirichlet = BoundaryConditions(
    BoundaryCondition("bottom", [0.0, 0.0, 0.0])
  )

  (; width, long, velocity₀, hrefinement, prefinement, t_end, Δt, dirichlet)
end

metrics, uh = solve_problem(problem_data)

Wdiff = metrics.Dvis + metrics.KE + metrics.EE - metrics.Wext

@show sum(metrics.Dvis) * problem_data.Δt
@show sum(Wdiff) * problem_data.Δt

let p = areaplot(metrics.time, [metrics.Dvis metrics.KE metrics.EE], label=["Dissipation" "Kinetic" "Elastic"])
  plot!(metrics.time, metrics.Wext, label="External", color=:black, lw=3, style=:dot)
  plot!(xlabel="Time [s]", ylabel="Energy [J]", legend=:outertop, legend_column=-1)
  display(p)
end
