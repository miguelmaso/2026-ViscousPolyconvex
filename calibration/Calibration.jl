using HyperFEM, HyperCalibration
using Optim
using CSV
using Plots


## Load data

quasi_static_set = CSV.read(joinpath(@__DIR__, "data", "quasi-static.csv"), UniaxialQuasiStaticTest)
cyclic_loading_set = CSV.read(joinpath(@__DIR__, "data", "cyclic_loading.csv"), UniaxialCyclicLoadingTest)


## Equilibrium component

build_model(μ) = NeoHookean3D(μ=μ, λ=0.0)
pn = [  "μ"]  # Parameter names
p0 = [  1e4]  # Initial seed

f(p) = loss(build_model, p, quasi_static_set)

eq_result = optimize(f, p0, NelderMead())
equilibrium_term = build_model(eq_result.minimizer...)
plot(equilibrium_term, quasi_static_set[1], xlabel="Stretch [-]", ylabel="Stress [KPa]", units_scale=1e-3)


## Non-equilibrium branches

build_branch(μ, τ) = ViscousPolyconvex(μ=μ, τ=τ)

