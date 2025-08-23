# SPACParams.jl

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jl-pkgs.github.io/SPACParams.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jl-pkgs.github.io/SPACParams.jl/dev)
[![CI](https://github.com/jl-pkgs/SPACParams.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/jl-pkgs/SPACParams.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/jl-pkgs/SPACParams.jl/branch/master/graph/badge.svg)](https://app.codecov.io/gh/jl-pkgs/SPACParams.jl/tree/master)

A Julia package for model parameter management and optimization.

## Features

- **Parameter Management**: Extract, update, and manage model parameters with metadata (bounds, units)
- **Parameter Optimization**: Optimize subsets of model parameters using the SCEUA algorithm
- **Goodness-of-Fit Functions**: Evaluate model performance with NSE, KGE, and other metrics
- **Flexible Model System**: Support for complex nested model structures

## Key Functions

- `Params(model)`: Extract all parameters from a model with their metadata
- `update!(model, names, values)`: Update specific model parameters by name
- `optimize_params!(model, objective_fn, param_names)`: Optimize selected parameters to minimize an objective function
- `GOF(obs, sim)`: Calculate goodness-of-fit metrics between observed and simulated data

## Parameter Optimization Example

```julia
using ModelParams

# Define your model with @bounds and @units macros
@bounds @units @with_kw mutable struct MyModel{FT} <: AbstractModel{FT}
  α::FT = 0.06 | (0.01, 0.10) | "units"
  β::FT = 2.0 | (1.0, 3.0) | "units"
end

model = MyModel{Float64}()

# Define objective function to minimize
objective_fn = m -> (m.α - 0.08)^2 + (m.β - 1.5)^2

# Optimize selected parameters
param_names = [:α, :β]
bounds = [(0.01, 0.10), (1.0, 3.0)]
optimal_params, cost, exitflag = optimize_params!(
  model, objective_fn, param_names;
  bounds=bounds, maxn=500
)
```

See [`examples/parameter_optimization.md`](examples/parameter_optimization.md) for a complete example.

> Dongdong Kong
