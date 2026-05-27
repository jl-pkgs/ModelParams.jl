# 检查某字段是否在父结构体上明确定义了 bound = nothing（即显式隐藏）
# ModelParamsdata 对有 `| val` 的字段生成特化方法 sig[2] == Type{<:T}，
# 无注解字段则退回到泛型方法 sig[2] == Type
function is_explicitly_hidden(T::Type, field)
    field isa Symbol || return false
    m = which(_meta, Tuple{Type{T},Val{field},Val{:bounds}})
    sig = Base.unwrap_unionall(m.sig)
    sig.parameters[2] !== Type && isnothing(_meta(T, Val{field}(), Val{:bounds}()))
end

# 把 bounds 分解成字段路径和对应的约束
function split_bounds(x::S) where {S}
    function categorize(field)
        is_explicitly_hidden(S, field) && return :skip
        value = getfield(x, field)
        T = typeof(value)
        # 不递归进 Array/Tuple 等内建类型，直接当普通字段处理（无 bound 则过滤掉）
        (T <: AbstractArray || T <: Tuple) && return :macro
        (has_definedbounds(value) || isstructtype(T)) ? :predef : :macro
    end
    predef, macro_fields = Symbol[], Symbol[]
    for f in fieldnames(S)
        cat = categorize(f)
        cat == :predef && push!(predef, f)
        cat == :macro  && push!(macro_fields, f)
        # :skip 的字段完全略过
    end
    predef, macro_fields
end


function get_params(x::S; path=[], with_unit=true) where {S}
    fs_predef, fs_macro = split_bounds(x)

    res_predef = map(field -> begin
            value = getfield(x, field)
            get_params(value; path=[path..., field], with_unit)
        end, fs_predef)

    res_macro = map(field -> begin
            value = getfield(x, field)
            unit = with_unit ? units(S, field) : ""
            (; path=[path..., field], name=field,
                value, type=eltype(value), bound=bounds(x, field), unit)
        end, fs_macro)
    res = vcat(res_macro, reduce(vcat, res_predef; init=eltype(res_macro)[]))
    filter(x -> !isnothing(x.bound), res)
end


function update!(model::S, paths::AbstractVector, values::Vector{FT},
    ; params::Union{Nothing,DataFrame}=nothing) where {S,FT}
    isnothing(params) && (params = parameters(model))

    path_idx = Dict(p => i for (i, p) in enumerate(params.path))
    length(path_idx) == nrow(params) || error("Duplicated parameters found in params!")
    @inbounds for (path, value) in zip(paths, values)
        idx = get(path_idx, path, nothing)
        isnothing(idx) && error("Parameter path $(path) not found in model!")
        update!(model, params.path[idx], value; type=params.type[idx])
    end
end

@inbounds function update!(model::S, path::AbstractVector, value::FT; type::Type) where {S,FT}
    if length(path) == 1
        setfield!(model, path[1], type(value))
    elseif length(path) > 1
        submodel = getfield(model, path[1])

        if isa(submodel, Vector) # 如果是多模型
            i = path[2]

            if typeof(submodel[i]) == FT
                submodel[i] = type(value)
                return
            end
            # 下面是应对Struct Vector
            update!(submodel[i], @view(path[3:end]), value; type)
        else
            update!(submodel, @view(path[2:end]), value; type)
        end
    end
end

function parameters(model; paths=nothing, with_unit=true)::DataFrame
    params = get_params(model; with_unit) |> DataFrame
    if !isnothing(paths)
        inds = indexin(paths, params.path)
        params = params[inds, :]
    end

    length(unique(params.unit)) == 1 && (params = params[:, Not(:unit)])
    return params
end

function get_opt_info(model; paths=nothing)
    df = parameters(model; paths)
    x0 = Float64.(df.value)
    lb = Float64[b[1] for b in df.bound]
    ub = Float64[b[2] for b in df.bound]
    _paths = df.path
    return x0, lb, ub, _paths
end

export get_params, parameters, update!
export get_opt_info
