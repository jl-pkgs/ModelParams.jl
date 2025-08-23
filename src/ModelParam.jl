import FieldMetadata: @bounds, bounds, @units, units


abstract type AbstractModel{FT} end


function unlist(list::Vector)
  res = []
  for x in list
    isa(x, Vector) ? append!(res, unlist(x)) : push!(res, x)
  end
  return map(x -> x, res)
end


_fieldname(x, field::Symbol) = field

"递归获取模型参数属性"
function get_ModelParamRecursive(x::T, path=Vector{Symbol}(); fun=bounds) where {FT,T<:AbstractModel{FT}}
  fileds = fieldnames(T) |> collect
  map(field -> begin
      value = getfield(x, field)
      _path = [path..., field]
      r = fun(x, field)
      typeof(value) == FT ? (_path => r) : get_ModelParamRecursive(value, _path; fun)
    end, fileds) |> unlist
end


function Params(model::AbstractModel)
  _names = get_ModelParamRecursive(model; fun=_fieldname)
  _bounds = get_ModelParamRecursive(model; fun=bounds)
  _values = get_ModelParamRecursive(model; fun=getfield)
  _units = get_ModelParamRecursive(model; fun=units)

  # 返回NamedTuple向量，符合Tables接口
  [(path=first(n), name=last(n), value=last(v),
    unit=last(u), bound=last(b))
   for (n, v, u, b) in zip(_names, _values, _units, _bounds)]
end


function update!(model::AbstractModel{FT}, names::Vector{Symbol}, values::Vector{FT},
  ; params::Union{Nothing,Vector{<:NamedTuple}}=nothing) where {FT}
  isnothing(params) && (params = Params(model))

  for (name, value) in zip(names, values)
    rows = filter(row -> row.name == name, params)
    @assert length(rows) == 1 "Duplicated parameters are not allowed!"
    update!(model, rows[1].path, value)
  end
end

function update!(model::AbstractModel{FT}, path::Vector{Symbol}, value::FT) where {FT}
  if length(path) == 1
    setfield!(model, path[1], value)
  elseif length(path) > 1
    submodel = getfield(model, path[1])
    update!(submodel, path[2:end], value)
  end
end


"""
    optimize_params!(model::AbstractModel, objective_fn::Function, param_names::Vector{Symbol}; 
                     bounds::Union{Nothing,Vector{Tuple{Real,Real}}}=nothing,
                     optimization_options...)

Optimize a subset of model parameters using the SCEUA algorithm.

# Arguments
- `model`: The model to optimize
- `objective_fn`: Function that takes the model and returns a cost value to minimize
- `param_names`: Vector of parameter names to optimize (must exist in model)
- `bounds`: Optional vector of (min, max) tuples for each parameter. If not provided, 
           attempts to extract from model metadata (may not work correctly in all cases)
- `optimization_options...`: Additional options passed to sceua() (maxn, kstop, etc.)

# Returns
- `(optimal_params, optimal_cost, exitflag)`: Optimized parameter values, final cost, and exit status

# Example
```julia
model = Photosynthesis_Rong2018{Float64}()
objective_fn = model -> sum(abs2, [model.α - 0.08, model.kQ - 0.5])  # minimize distance to target
param_names = [:α, :kQ]
bounds = [(0.01, 0.10), (0.10, 1.0)]  # Optional manual bounds
optimal_params, cost, exitflag = optimize_params!(model, objective_fn, param_names; bounds, maxn=1000)
```
"""
function optimize_params!(model::AbstractModel{FT}, objective_fn::Function, param_names::Vector{Symbol}; 
                         bounds::Union{Nothing,Vector{<:Tuple{<:Real,<:Real}}}=nothing,
                         kw...) where {FT}
  # Extract all model parameters
  params = Params(model)
  
  # Filter to get only the parameters we want to optimize
  opt_params = [p for p in params if p.name in param_names]
  
  # Check that all requested parameters exist
  found_names = [p.name for p in opt_params]
  missing_names = setdiff(param_names, found_names)
  if !isempty(missing_names)
    error("Parameters not found in model: $(missing_names)")
  end
  
  # Check for duplicates
  if length(unique(found_names)) != length(found_names)
    error("Duplicate parameter names found: $(found_names)")
  end
  
  # Sort opt_params to match the order of param_names
  param_order = Dict(name => i for (i, name) in enumerate(param_names))
  sort!(opt_params, by = p -> param_order[p.name])
  
  # Extract bounds and initial values
  n_params = length(opt_params)
  x0 = FT.(getfield.(opt_params, :value))
  
  if bounds !== nothing
    # Use provided bounds
    if length(bounds) != n_params
      error("Number of bounds ($(length(bounds))) must match number of parameters ($n_params)")
    end
    bl = FT.([b[1] for b in bounds])
    bu = FT.([b[2] for b in bounds])
  else
    # Try to extract bounds from model metadata
    bl = FT.([p.bound[1] for p in opt_params])
    bu = FT.([p.bound[2] for p in opt_params])
    
    # Check if bounds extraction failed (all (0,1) suggests it failed)
    if all(b == (0.0, 1.0) for b in zip(bl, bu))
      @warn "Bounds extraction may have failed - all bounds are (0.0, 1.0). Consider providing manual bounds."
    end
  end
  
  # Create wrapper function for SCEUA
  function sceua_objective(x::Vector{FT})
    # Update model with new parameter values
    update!(model, param_names, x; params)
    # Evaluate objective function
    cost = objective_fn(model)
    return convert(FT, cost)
  end
  
  # Run optimization
  optimal_x, optimal_cost, exitflag = sceua(sceua_objective, x0, bl, bu; kw...)
  
  # Update model with optimal parameters
  update!(model, param_names, optimal_x; params)
  
  # Return results as named tuple for clarity
  optimal_params = [(name=param_names[i], value=optimal_x[i]) for i in 1:n_params]
  
  return optimal_params, optimal_cost, exitflag
end
