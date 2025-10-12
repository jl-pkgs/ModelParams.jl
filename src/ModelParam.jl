import FieldMetadata: @metadata, @units, units

@metadata bounds nothing

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
function get_ModelParamRecur(x::T, path=Vector{Symbol}(); fun=bounds) where {FT,T<:AbstractModel{FT}}
  fileds = fieldnames(T) |> collect
  map(field -> begin
      value = getfield(x, field)
      _path = [path..., field]
      r = fun(x, field)
      TYPE = typeof(value)

      if TYPE == FT
        return (_path => r)
      elseif TYPE <: AbstractModel
        get_ModelParamRecur(value, _path; fun)
      elseif isMultiModels(value)
        map(i -> begin
            _path2 = [_path..., i]
            get_ModelParamRecur(value[i], _path2; fun)
        end, 1:length(value))
      else
        return []
      end
    end, fileds) |> unlist
end


isMultiModels(value) = isa(value, Vector) && eltype(typeof(value)) <: AbstractModel

function Params(model::AbstractModel; na_rm::Bool=true)
  _names = get_ModelParamRecur(model; fun=_fieldname)
  _bounds = get_ModelParamRecur(model; fun=bounds)
  _values = get_ModelParamRecur(model; fun=getfield)
  _units = get_ModelParamRecur(model; fun=units)

  # 返回NamedTuple向量，符合Tables接口
  params = [
    (name=last(n), value=last(v), bound=last(b), unit=last(u), path=first(n))
    for (n, v, u, b) in zip(_names, _values, _units, _bounds)] |> DataFrame

  lgl = (map(x -> !isnothing(x), params.bound)) # bounds不为NaN的params
  na_rm && (params = params[lgl, :])
  return params
end


function update!(model::AbstractModel{FT}, paths::Vector, values::Vector{FT},
  ; params::Union{Nothing,DataFrame}=nothing) where {FT}
  isnothing(params) && (params = Params(model))

  for (path, value) in zip(paths, values)
    rows = filter(row -> row.path == path, params)
    @assert size(rows, 1) == 1 "Duplicated parameters are not allowed!"
    update!(model, rows.path[1], value)
  end
end

function update!(model::AbstractModel{FT}, path::Vector, value::FT) where {FT}
  if length(path) == 1
    setfield!(model, path[1], value)
  elseif length(path) > 1
    submodel = getfield(model, path[1]) # 
    if isMultiModels(submodel) # 如果是多模型
      models = submodel
      i = path[2]
      update!(models[i], path[3:end], value)
    else
      update!(submodel, path[2:end], value)
    end
  end
end

function get_bound(bound::Vector)
  lower = map(x -> x[1], bound)
  upper = map(x -> x[2], bound)
  lower, upper
end

get_bound(params::DataFrame) = get_bound(params.bound)
