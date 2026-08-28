using Revise
using HyperFEM, HyperCalibration
using Optim
using CSV
using Plots

default(
  linewidth = 2,
  mswidth = 0,
  palette = :seaborn_colorblind
)

## Load data

quasi_static_set = CSV.read(joinpath(@__DIR__, "data", "quasi-static.csv"), UniaxialQuasiStaticTest)
cyclic_loading_set = CSV.read(joinpath(@__DIR__, "data", "cyclic-loading.csv"), UniaxialCyclicLoadingTest)

sort!(cyclic_loading_set, by = r -> rate(r))
sort!(cyclic_loading_set, by = r -> max_stretch(r))


## Equilibrium component

build_equilibrium(μ) = NeoHookean3D(μ=μ, λ=0.0)
pn = [  "μ"]  # Parameter names
p0 = [  1e4]  # Initial seed

f(p) = normalized_mse(build_equilibrium, p, quasi_static_set)

equil_opt = optimize(f, p0, NelderMead())
equil_result = CalibrationResult(build_equilibrium, Optim.minimizer(equil_opt), quasi_static_set)

println(parameter_stats(equil_result, names=pn))
display(plot(equil_result, quasi_static_set[1], xlabel="Stretch [-]", ylabel="Stress [KPa]", units_scale=1e-3))


## Non-equilibrium branches

build_branch(μ, τ) = ViscousPolyconvex(μ=μ, τ=τ)
build_branches(p...) = map(splat(build_branch), Iterators.partition(p, 2))
build_visco(p...) = GeneralizedMaxwell(equil_result.model, build_branches(p...)...)
n_branches = 2
pn = reduce(vcat, ["μ$i", "τ$i"] for i in 1:n_branches)  # Parameter names
p0 = reduce(vcat, [  1e4,   1.0] for _ in 1:n_branches)  # Initial seed
lb = reduce(vcat, [  1e3,  -1.0] for _ in 1:n_branches)  # Lower search limits
ub = reduce(vcat, [  1e5,   4.0] for _ in 1:n_branches)  # Upper search limits

f(p) = normalized_mse(build_visco, p, cyclic_loading_set)

noneq_opt = optimize(f, p0, ParticleSwarm())
noneq_result = CalibrationResult(build_visco, Optim.minimizer(noneq_opt), cyclic_loading_set)

display(MIME("text/latex"), parameter_stats(noneq_result, names=pn))
fixed_stretch_subset = filter(r -> max_stretch(r) ≈ 1.98, cyclic_loading_set)
fixed_loading_rate_subset = filter(r -> rate(r) ≈ 0.03, cyclic_loading_set)
display(plot(noneq_result, fixed_stretch_subset, xlabel="Stretch [-]", ylabel="Stress [KPa]", units_scale=1e-3, labels=map(pretty_label(rate), fixed_stretch_subset)))
display(plot(noneq_result, fixed_loading_rate_subset, xlabel="Stretch [-]", ylabel="Stress [KPa]", units_scale=1e-3, labels=map(pretty_label(max_stretch), fixed_loading_rate_subset)))


## Uncertainty visualization

rand_params = sample_parameters(noneq_result)
models = map(splat(build_visco), eachcol(rand_params))

experiment = first(filter(r -> max_stretch(r) ≈ 1.98 && rate(r) ≈ 0.03, cyclic_loading_set))  # 1.98
p = plot(xlabel="Stretch [-]", ylabel="Stress [KPa]")
plot!(models, experiment, color=1, alpha=0.05, label=false, units_scale=1e-3)
plot!(noneq_result, experiment, color=[1 :black], label=["Model" "Data"], units_scale=1e-3)
display(p);

