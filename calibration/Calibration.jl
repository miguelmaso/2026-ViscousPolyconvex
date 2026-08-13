using HyperFEM
using HyperCalibration
using CSV
using Plots

experiments = CSV.read(joinpath(@__DIR__, "data", "quasi-static.csv"), UniaxialQuasiStaticTest)

build_model(μ, N) = NeoHookean3D(μ=μ, λ=0.0)
pn = [  "μ"]  # Parameter names
p0 = [  1e4]  # Initial seed

f(p) = loss(build_model, p, experiments)

result = optimize(f, p0, NelderMead())
model = build_model(result.minimizer...)
plot(model, experiments[1], xlabel="Stretch [-]", ylabel="Stress [KPa]", units_scale=1e-3)
