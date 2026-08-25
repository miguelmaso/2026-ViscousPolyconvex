# 2026-ViscousPolyconvex

Numerical examples showcasing the use of the viscous model implemented in [HyperFEM](https://github.com/MultiSimOLab/HyperFEM.jl) ([viscous polyconvex](https://github.com/MultiSimOLab/HyperFEM.jl/blob/main/src/PhysicalModels/ViscousPolyconvex.jl)).

## VHB 4905 characterization

Parameters fitted with the [HyperCalibration](https://github.com/miguelmaso/HyperCalibration.jl) library:

Parameter | Value
--------- | -----
μe        | 13.4e3  Pa
κr        | 2.5e6   Pa
μ1        | 30.9e3  Pa
τ1        | 6.65    s
μ2        | 11.4e3  Pa
τ2        | 144.1   s

![](docs/img/cyclic_loading_fixed_rate.png)


## Numerical example

Dynamic integration of a thin membrane with out-of-plane displacement.

![](docs/img/ViscoTimeIntegrator.gif)
