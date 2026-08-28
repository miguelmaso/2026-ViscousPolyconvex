# 2026-ViscousPolyconvex

Numerical examples showcasing the use of the viscous model implemented in [HyperFEM](https://github.com/MultiSimOLab/HyperFEM.jl) ([viscous polyconvex](https://github.com/MultiSimOLab/HyperFEM.jl/blob/main/src/PhysicalModels/ViscousPolyconvex.jl)).

## VHB 4905 characterization

The material characterization has been performed with the [HyperCalibration](https://github.com/miguelmaso/HyperCalibration.jl) library, and the data is described in the [data](https://github.com/miguelmaso/2026-ViscousPolyconvex/tree/main/calibration/data) folder. A new-Hookean constitutive model has been selected both for the equilibrium and the non-equilibrium terms, and the fitted parameters are summarized below:

Parameter | Estimate ± Margin  | Rel. Err (%) | Sensitivity
----------|--------------------|--------------|----------
μe        | 1.34e+04 ± 2.4e+02 | 1.8          | 1278.6
μ1        | 3.12e+04 ± 3.4e+03 | 11.0         | 18.0
t1        |     6.45 ± 1.9     | 29.0         | 12.1
μ2        | 1.16e+04 ± 3.2e+03 | 27.5         | 48.1
t2        |      140 ± 76      | 54.2         | 6.0

![](docs/img/cyclic_loading_fixed_rate.png)


## Numerical example

Dynamic integration of a thin membrane with out-of-plane displacement. The previously characterized constitutive model is used for the simulation.

![](docs/img/ViscoTimeIntegrator.gif)
