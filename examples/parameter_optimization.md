# Parameter Optimization Example

This example demonstrates how to use the `optimize_params!` function to optimize a subset of model parameters.

```julia
using ModelParams, Parameters

# Define model structures (same as in the tests)
abstract type AbstractWaterConsGPPModel{FT} <: AbstractModel{FT} end
abstract type AbstractPhotosynthesisModel{FT} <: AbstractModel{FT} end

@bounds @units @with_kw mutable struct β_GPP_Zhang2019{FT} <: AbstractWaterConsGPPModel{FT}
  VPDmin::FT = 0.9 | (0.65, 1.5) | "kPa"
  VPDmax::FT = 4.0 | (3.50, 6.5) | "kPa"
end

@bounds @units @with_kw mutable struct Photosynthesis_Rong2018{FT} <: AbstractPhotosynthesisModel{FT}
  α::FT = 0.06 | (0.01, 0.10) | "μmol CO2 [μmol PAR]⁻¹"
  η::FT = 0.04 | (0.01, 0.07) | "μmol m⁻² s⁻¹ [μmol m⁻² s⁻¹]⁻¹"
  VCmax25::FT = 50.00 | (5.00, 120.00) | "μmol m⁻² s⁻¹"
  d_pc::FT = 2.0 | (0.0, 5.0) | "-"
  kQ::FT = 0.45 | (0.10, 1.0) | "-"
  watercons::AbstractWaterConsGPPModel{FT} = β_GPP_Zhang2019{FT}()
end

# Create a model instance
FT = Float64
model = Photosynthesis_Rong2018{FT}()

# Display all model parameters
println("All model parameters:")
params = Params(model)
for p in params
  println("  $(p.name): $(p.value)")
end

# Define an objective function to minimize
# Example: minimize squared distance from target values
function objective_function(m)
  target_α = 0.08
  target_kQ = 0.5
  target_VCmax25 = 60.0
  
  return (m.α - target_α)^2 + 
         (m.kQ - target_kQ)^2 + 
         (m.VCmax25 - target_VCmax25)^2
end

println("\nInitial objective value: $(objective_function(model))")

# Optimize only a subset of parameters
param_names = [:α, :kQ, :VCmax25]
manual_bounds = [
  (0.01, 0.10),    # bounds for α
  (0.10, 1.0),     # bounds for kQ  
  (5.00, 120.00)   # bounds for VCmax25
]

println("\nOptimizing parameters: $param_names")

# Run optimization
optimal_params, optimal_cost, exitflag = optimize_params!(
  model, objective_function, param_names;
  bounds=manual_bounds,
  maxn=500,
  verbose=false
)

println("\nOptimization results:")
println("  Exit flag: $exitflag")
println("  Final cost: $optimal_cost")
println("  Optimized parameters:")
for p in optimal_params
  println("    $(p.name): $(p.value)")
end

println("\nFinal model state:")
for p in params
  println("  $(p.name): $(getfield(model, p.name))")
end
```

## Key Features

1. **Selective Parameter Optimization**: Choose which parameters to optimize using `param_names`
2. **Manual Bounds**: Provide custom bounds with the `bounds` parameter (recommended)
3. **Automatic Bounds**: The function attempts to extract bounds from model metadata, but this may not work correctly in all cases
4. **SCEUA Integration**: Uses the robust SCEUA optimization algorithm
5. **Model Update**: The model is automatically updated with the optimal parameter values

## Usage Notes

- It's recommended to provide manual bounds for better control and reliability
- The function preserves the order of parameters as specified in `param_names`
- All SCEUA optimization options (maxn, kstop, pcento, etc.) can be passed as keyword arguments
- The function validates that all requested parameters exist in the model