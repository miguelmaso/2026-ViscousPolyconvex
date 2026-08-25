using Revise
using HyperFEM, HyperFEM.ComputationalModels.PostMetrics
using Gridap, Gridap.FESpaces, Gridap.Geometry
using GridapSolvers, GridapSolvers.NonlinearSolvers
using Printf
using JLD2


## Domain

function generate_tessellation(; width, thick, hrefinement, hdivisions, args...)
  domain = (-0.5width, 0.5width, -0.5width, 0.5width, 0.0, thick)
  partition = hrefinement .* (hdivisions, hdivisions, 1)
  geometry = CartesianDiscreteModel(domain, partition)
  labels = get_face_labeling(geometry)
  add_tag_from_tags!(labels, "top",    CartesianTags.faceXY1⁺)
  add_tag_from_tags!(labels, "bottom", CartesianTags.faceXY0⁺)
  add_tag_from_vertex_filter!(labels, "center", geometry, p -> abs(p[1]) < 0.51width/hdivisions && abs(p[2]) < 0.51width/hdivisions)
  geometry
end


## Constitutive model

function build_model(; args...)

  # Equilibrium branch
  μe = 13.4e3  # [Pa]
  κr = 2.5e6   # [Pa]

  # Non-equilibrium branches
  μ1 = 30.9e3  # [Pa]
  τ1 = 6.65    # [s]
  μ2 = 11.4e3  # [Pa]
  τ2 = 144.1   # [s]

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

  writevtk(geometry, outpath)

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
  υh  = FEFunction(Uυ,  zero_free_values(Uυ))

  Fh  = F∘∇(uh⁺)'
  Fh⁻ = F∘∇(uh⁻)'
  Ah  = CellState(model, dΩ)
  Ph⁻ = CellState(zero(TensorValue{3,3,Float64}), dΩ)

  # Weak forms: residual and jacobian

  update_time_step!(model, Δt)
  Ψ, ∂Ψ∂F, ∂∂Ψ∂FF = model()
  D = Dissipation(model)
  P_updater(P, F, Fn, A...) = (true, ∂Ψ∂F(F, Fn, A...))

  # res_mec(time) = (u, v) -> ∫( ∇(v)' ⊙ (∂Ψ∂F ∘ (F∘(∇(u)'), Fh⁻, Ah...)) )dΩ
  # jac_mec(time) = (u, du, v) -> ∫( ∇(v)' ⊙ ((∂∂Ψ∂FF ∘ (F∘(∇(u)'), Fh⁻, Ah...)) ⊙ ∇(du)') )dΩ

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

  fields = (:time, :Ψmec, :Ψdir, :Dvis)
  metrics = NamedTuple(f => Float64[] for f in fields)

  function post_metrics!(data, step, time)
    b = assemble_vector(v -> residual(uh⁺, v), DirichletFESpace(Vu))[:]
    dudt_dir = (get_dirichlet_dof_values(Uu) - get_dirichlet_dof_values(Uu⁻)) / Δt
    push!(data.time, time)
    push!(data.Ψmec, sum(residual(uh⁺, uh⁺-uh⁻))/Δt)
    push!(data.Ψdir, b · dudt_dir)
    push!(data.Dvis, sum(∫( D∘(Fh, Fh⁻, Ah...) )dΩ))
  end

  function post_vtk!(pvd, step, time)
    if mod(step, 10) == 0
      Ph = interpolate_L2_field(∂Ψ∂F ∘ (Fh, Fh⁻, Ah...), Ω, dΩ)
      Jh = interpolate_L2_field(J∘Fh, Ω, dΩ)
      pvd[time] = createvtk(Ω, outpath * @sprintf("_%03d", step), cellfields=["u" => uh⁺, "J" => Jh, "P" => Ph])
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
        printstyled(@sprintf("Step: %i\nTime: %.3f s\n", step, time), color=:green, bold=true)

        #-----------------------------------------
        # Apply boundary conditions and solve
        #-----------------------------------------
        TrialFESpace!(Uu, data.dirichlet_u, time)

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
  width = 0.2
  thick = 0.001
  speed = 0.1
  hdivisions = 35
  hrefinement = 1
  prefinement = 2
  t_end = 0.5
  CFL = 0.2
  Δt = CFL * thick / (prefinement * hrefinement * speed)

  dir_u_tags = ["center"]
  dir_u_values = [[0.0, 0.0, 1.0]]
  dir_u_time = [t -> t*speed]
  dirichlet_u = DirichletBC(dir_u_tags, dir_u_values, dir_u_time)

  (; width, thick, speed, hrefinement, hdivisions, prefinement, t_end, Δt, dirichlet_u)
end

solve_problem(problem_data)
