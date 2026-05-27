export parameters, update!
export AbstractLayers, MultiLayer, Layers

using DataFrames

# 抽象基类：FT 元素类型，N 层数（去掉旧设计中的 S 参数）
abstract type AbstractLayers{FT,N} end
# abstract type AbstractModel{FT} end

# 统一的 SoA 多层容器
# - S:  标量结构体类型（如 Campbell{FT}）
# - NT: 底层 NamedTuple 类型（由 S 的字段集和 FT 决定，构造时自动推断）
struct MultiLayer{FT,N,S,NT<:NamedTuple} <: AbstractLayers{FT,N}
    data::NT
end

# 把 kwargs 值规范化为 Vector{FT}
_ml_field(::Type{FT}, N::Int, v::AbstractVector) where {FT} = convert(Vector{FT}, v)
_ml_field(::Type{FT}, N::Int, v::Number) where {FT} = fill(FT(v), N)

# @generated 构造器：从 S 的字段集自动展开 SoA NamedTuple
@generated function MultiLayer{FT,N,S}(; kwargs...) where {FT,N,S}
    fnames = fieldnames(S)
    ftypes = fieldtypes(S)
    keep = Symbol[fnames[i] for (i, t) in enumerate(ftypes) if t <: AbstractFloat]
    keep_tuple = Tuple(keep)

    s_default = S()
    defaults = Float64[Float64(getfield(s_default, n)) for n in keep]

    val_exprs = map(zip(keep, defaults)) do (n, d)
        qn = QuoteNode(n)
        :(haskey(kwargs, $qn) ? _ml_field($FT, $N, kwargs[$qn]) : fill($FT($d), $N))
    end

    quote
        data = NamedTuple{$keep_tuple}(tuple($(val_exprs...)))
        MultiLayer{$FT,$N,$S,typeof(data)}(data)
    end
end

# 字段访问：转发到底层 NamedTuple
@inline function Base.getproperty(x::MultiLayer, name::Symbol)
    name === :data && return getfield(x, :data)
    return getproperty(getfield(x, :data), name)
end

Base.propertynames(::MultiLayer{FT,N,S}) where {FT,N,S} = (:data, fieldnames(S)...)

# AoS 单层：根据 NT 的字段名重建标量 S
@generated function Base.getindex(x::MultiLayer{FT,N,S,NT}, i::Int) where {FT,N,S,NT}
    fnames = NT.parameters[1]
    kw_exprs = [:($n = x.data.$n[i]) for n in fnames]
    :(S(; $(kw_exprs...)))
end

Base.length(::MultiLayer{FT,N}) where {FT,N} = N

function Base.Vector(x::MultiLayer{FT,N,S}) where {FT,N,S}
    return S[x[i] for i in 1:N]
end

# 从单层 AoS 实例构造 SoA MultiLayer
@generated function Layers(p::P, N::Int) where {P}
    FT = Float64
    fnames = fieldnames(P)
    ftypes = fieldtypes(P)
    for ft in ftypes
        ft <: AbstractFloat && (FT = ft; break)
    end
    float_fields = [fnames[i] for (i, t) in enumerate(ftypes) if t <: AbstractFloat]
    kwargs = [:($f = fill(getfield(p, $(QuoteNode(f))), N)) for f in float_fields]
    quote
        MultiLayer{$FT,N,$P}(; $(kwargs...))
    end
end


has_definedbounds(x) = false
has_definedbounds(x::AbstractLayers) = true

# 专属 update!：MultiLayer 的字段（如 :θ_sat）不是真实字段，需走 getproperty
function update!(x::MultiLayer, path::Vector, value::FT; type::Type) where {FT}
    if length(path) >= 2
        field = path[1]
        idx = path[2]
        vec = getproperty(x, field)
        vec[idx] = type(value)
        return
    end
    error("update! on MultiLayer requires path of length >= 2, got $path")
end

# 多层参数收集：从 S 取字段集和 bounds/units 元数据
@generated function _ml_float_fields(::Type{S}) where {S}
    keep = [f for (f, t) in zip(fieldnames(S), fieldtypes(S)) if t <: AbstractFloat]
    :($(Tuple(keep)))
end

function get_params(x::MultiLayer{FT,N,S}; path=[], with_unit=true) where {FT,N,S}
    fnames = _ml_float_fields(S)
    isempty(fnames) && return []
    res = map(fnames) do field
        value = getproperty(x, field)
        _path = [path..., field]
        bound = bounds(S, field)
        unit = with_unit ? units(S, field) : ""
        [(; path=[_path..., i], name=field,
            value=value[i], type=eltype(value), bound, unit) for i in 1:N]
    end
    filter(r -> !isnothing(r.bound), reduce(vcat, res))
end


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
    fields = fieldnames(S)
    (
        filter(f -> categorize(f) == :predef, fields),
        filter(f -> categorize(f) == :macro, fields)
        # :skip 的字段完全略过
    )
end


function get_params(x::S; path=[], with_unit=true) where {S}
    fs_predef, fs_macro = split_bounds(x)

    res_predef = map(field -> begin
            value = getfield(x, field)
            get_params(value; path=[path..., field])
        end, fs_predef)

    res_macro = map(field -> begin
            # @show bounds(x, field)
            value = getfield(x, field)
            unit = with_unit ? units(S, field) : ""
            (; path=[path..., field], name=field,
                value, type=eltype(value), bound=bounds(x, field), unit)
        end, fs_macro)
    res = vcat(res_macro..., res_predef...)
    filter(x -> !isnothing(x.bound), res)
end


function update!(model::S, paths::Vector, values::Vector{FT},
    ; params::Union{Nothing,DataFrame}=nothing) where {S,FT}
    isnothing(params) && (params = parameters(model))

    path_idx = Dict(p => i for (i, p) in enumerate(params.path))
    length(path_idx) == nrow(params) || error("Duplicated parameters found in params!")
    for (path, value) in zip(paths, values)
        idx = get(path_idx, path, nothing)
        isnothing(idx) && error("Parameter path $(path) not found in model!")
        update!(model, params.path[idx], value; type=params.type[idx])
    end
end

function update!(model::S, path::Vector, value::FT; type::Type) where {S,FT}
    if length(path) == 1
        # @show model, path[1], value
        setfield!(model, path[1], type(value))
    elseif length(path) > 1
        submodel = getfield(model, path[1]) #

        if isa(submodel, Vector) # 如果是多模型
            models = submodel
            i = path[2]

            if typeof(models[i]) == FT
                models[i] = type(value)
                return
            end
            # 下面是应对Struct Vector
            update!(models[i], path[3:end], value; type)
        else
            update!(submodel, path[2:end], value; type)
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
    paths = df.path
    return x0, lb, ub, paths
end

export get_params, parameters, update!
