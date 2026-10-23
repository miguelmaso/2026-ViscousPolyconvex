using Revise
using HyperFEM, HyperCalibration
using Optim
using CSV
using Plots
using Random
using Printf
using LaTeXStrings

default(
  linewidth = 2,
  mswidth = 0,
  palette = :seaborn_colorblind,
  size = (390, 290),
  guidefontsize = 9,
  tickfontsize = 8,
  legendfontsize = 8,
  fontfamily="Computer Modern",
)
Random.seed!(1234)


## Configuration: EPU40, Hossain, Navaratne & Perić (2020), see data/EPU40/README.md

const MATERIAL    = "EPU40"
const DATA_DIR    = joinpath(@__DIR__, "data", MATERIAL)
const FIG_DIR     = mkpath(joinpath(@__DIR__, "..", "latex", "figures", MATERIAL))
const UNITS       = "MPa"
const UNITS_SCALE = 1e-6

# Equilibrium component. The multi-step relaxation data (id 2) reach λ = 5 and are
# concave in nominal stress, which a neo-Hookean term cannot follow.
build_equilibrium(μ1, μ2) = MooneyRivlin3D(μ1=μ1, μ2=μ2, λ=0.0)
const EQUIL_NAMES = ["μ1", "μ2"]
const EQUIL_SEED  = [5e5, 5e5]

# Non-equilibrium branches (μ, τ)
const N_BRANCHES  = 2
const BRANCH_SEED = [5e5, 1.0]

# Plots
const QS_PLOT        = [1 => "Quasi-static", 2 => "Multi-step relaxation"]  # id => label
const FIXED_STRETCH  = 3.0     # λmax of the fixed-stretch subset
const FIXED_RATE     = 0.128   # rate of the fixed-rate subset
const EXAMPLE        = (λmax = 3.0, rate = 0.128)  # test of the uncertainty plot
const EXAMPLE_YLIMS  = :auto
rate_label(d)    = @sprintf("%.3f /s", rate(d))
stretch_label(d) = pretty_label(max_stretch, d)
marker_stride(d) = max(1, round(Int, (max_stretch(d) - 1) / 2))  # ~85 markers per test, as λmax = 3

const XLABEL = "Stretch λ [-]"
const YLABEL = "Nominal stress $(L"P_{11}") [$(UNITS)]"


## Plot helper: model as a line, data as markers (every `stride` points, peak and last point kept)

function plot_fit!(p, model, data; color, marker_color=color, label="", data_label=false, stride=marker_stride(data))
  x = independent_variable(data)
  y_true, y_pred = experiment_prediction(model, data) .* UNITS_SCALE
  i = unique(sort([1:stride:length(x); argmax(x); length(x)]))
  plot!(p, x, y_pred; color, label)
  scatter!(p, x[i], y_true[i]; color=marker_color, label=data_label)
  p
end


## Load data

quasi_static_set = CSV.read(joinpath(DATA_DIR, "quasi-static.csv"), UniaxialQuasiStaticTest)
cyclic_loading_set = CSV.read(joinpath(DATA_DIR, "cyclic-loading.csv"), UniaxialCyclicLoadingTest)

sort!(quasi_static_set, by = r -> r.id)
sort!(cyclic_loading_set, by = r -> rate(r))
sort!(cyclic_loading_set, by = r -> max_stretch(r))


## Equilibrium component

f(p) = normalized_mse(build_equilibrium, p, quasi_static_set)

equil_opt = optimize(f, EQUIL_SEED, NelderMead())
equil_result = CalibrationResult(build_equilibrium, Optim.minimizer(equil_opt), quasi_static_set)

display(parameter_stats(equil_result, names=EQUIL_NAMES))
let p = plot(xlabel=XLABEL, ylabel=YLABEL)
  for (i, (id, label)) in enumerate(QS_PLOT)
    plot_fit!(p, equil_result.model, only(filter(r -> r.id == id, quasi_static_set)), color=i, label=label, stride=1)
  end
  display(p)
  savefig(p, joinpath(FIG_DIR, "calibration_equilibrium.pdf"))
end


## Non-equilibrium branches

build_branch(μ, τ) = ViscousPolyconvex(μ=μ, τ=τ)
build_branches(p...) = map(splat(build_branch), Iterators.partition(p, 2))
build_visco(p...) = GeneralizedMaxwell(equil_result.model, build_branches(p...)...)
pn = reduce(vcat, ["μ$i", "τ$i"] for i in 1:N_BRANCHES)  # Parameter names
p0 = repeat(BRANCH_SEED, N_BRANCHES)                      # Initial seed

f(p) = normalized_mse(build_visco, p, cyclic_loading_set)

noneq_opt = optimize(f, p0, ParticleSwarm())
noneq_result = CalibrationResult(build_visco, Optim.minimizer(noneq_opt), cyclic_loading_set)

display(parameter_stats(noneq_result, names=pn))
fixed_stretch_subset = filter(r -> max_stretch(r) ≈ FIXED_STRETCH, cyclic_loading_set)
fixed_loading_rate_subset = filter(r -> rate(r) ≈ FIXED_RATE, cyclic_loading_set)
let p = plot(xlabel=XLABEL, ylabel=YLABEL)
  for (i, d) in enumerate(fixed_stretch_subset)
    plot_fit!(p, noneq_result.model, d, color=i, label=rate_label(d))
  end
  display(p)
  savefig(p, joinpath(FIG_DIR, "calibration_fixed_stretch.pdf"))
end
let p = plot(xlabel=XLABEL, ylabel=YLABEL)
  for (i, d) in enumerate(fixed_loading_rate_subset)
    plot_fit!(p, noneq_result.model, d, color=i, label=stretch_label(d))
  end
  display(p)
  savefig(p, joinpath(FIG_DIR, "calibration_fixed_rate.pdf"))
end


## Uncertainty visualization

rand_params = sample_parameters(noneq_result)
models = map(splat(build_visco), eachcol(rand_params))

experiment = only(filter(r -> max_stretch(r) ≈ EXAMPLE.λmax && rate(r) ≈ EXAMPLE.rate, cyclic_loading_set))
let p = plot(xlabel=XLABEL, ylabel=YLABEL, ylims=EXAMPLE_YLIMS)
  plot!(models, experiment, color=1, alpha=0.05, label=false, units_scale=UNITS_SCALE)
  plot_fit!(p, noneq_result.model, experiment, color=1, marker_color=:black, label="Model", data_label="Data")
  display(p)
  savefig(p, joinpath(FIG_DIR, "uncertainty_$(N_BRANCHES)_branches.pdf"))
end


## Time-step check: the internal variables are integrated with Δt = Δλ/rate of the data.
## Re-evaluate the calibrated model with Δλ/4 and compare at the data points.

function refined(d::UniaxialCyclicLoadingTest, n=4)
  λ = stretches(d)
  λf = reduce(vcat, [range(λ[k-1], λ[k], length=n+1)[2:end] for k in 2:length(λ)]; init=[λ[1]])
  UniaxialCyclicLoadingTest(d.id, λf, rate(d), zeros(length(λf)))
end

for d in cyclic_loading_set
  _, y  = experiment_prediction(noneq_result.model, d)
  _, yf = experiment_prediction(noneq_result.model, refined(d))
  @printf("Time-step check  %s, %s: Δt = %.3f s, max change with Δλ/4 = %.2f %% of peak\n",
          stretch_label(d), rate_label(d), time_step(d), 100 * maximum(abs.(yf[1:4:end] .- y)) / maximum(abs.(y)))
end
