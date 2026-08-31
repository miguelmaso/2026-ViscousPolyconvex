
Library requirements:

- Neumann boundary conditions: HyperFEM.jl v0.0.6
- Dirichlet boundary condition: HyperFEM.jl v0.0.6 pre-release (see the sources section in Project.toml)

**Option A**: Pointing to local copy of HyperFEM
[sources]
HyperFEM = {path = "../../HyperFEM.jl"}

**Option B**: Pointing to github repo
[sources]
HyperFEM = {url = "[../../HyperFEM.jl](https://github.com/MultiSimOLab/HyperFEM.jl)"}
