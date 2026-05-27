export parameters, update!
export AbstractLayers, MultiLayer, Layers

using DataFrames

# 抽象基类：FT 元素类型，N 层数（去掉旧设计中的 S 参数）
abstract type AbstractLayers{FT,N} end

# 统一的 SoA 多层容器
# - S:  标量结构体类型（如 Campbell{FT}）
# - NT: 底层 NamedTuple 类型（由 S 的字段集和 FT 决定，构造时自动推断）
struct MultiLayer{FT,N,S,NT<:NamedTuple} <: AbstractLayers{FT,N}
    data::NT
end

# 把 kwargs 值规范化为 Vector{FT}
_ml_field(::Type{FT}, N::Int, v::AbstractVector) where {FT} = convert(Vector{FT}, v)
_ml_field(::Type{FT}, N::Int, v::Number) where {FT} = fill(FT(v), N)

# !注意：MultiLayer自动只保留了Float类型的字段
# @generated 构造器：从 S 的字段集自动展开 SoA NamedTuple
@generated function MultiLayer{FT,N,S}(; kwargs...) where {FT,N,S}
    fields = fieldnames(S)
    ftypes = fieldtypes(S)
    keep = Symbol[fields[i] for (i, t) in enumerate(ftypes) if t <: AbstractFloat]
    keep_tuple = Tuple(keep)

    s_default = S() # S must have a zero-argument constructor (e.g. via @with_kw)
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

Base.propertynames(::MultiLayer{FT,N,S}) where {FT,N,S} = _ml_float_fields(S)

# AoS 单层：根据 NT 的字段名重建标量 S
@generated function Base.getindex(x::MultiLayer{FT,N,S,NT}, i::Int) where {FT,N,S,NT}
    fields = fieldnames(NT)
    kw_exprs = [:($n = x.data.$n[i]) for n in fields]
    :(S(; $(kw_exprs...)))
end

Base.length(::MultiLayer{FT,N}) where {FT,N} = N

function Base.Vector(x::MultiLayer{FT,N,S}) where {FT,N,S}
    return S[x[i] for i in 1:N]
end

# 从单层 AoS 实例构造 SoA MultiLayer
@generated function Layers(p::P, N::Int) where {P}
    FT = Float64
    fields = fieldnames(P)
    ftypes = fieldtypes(P)
    for ft in ftypes
        ft <: AbstractFloat && (FT = ft; break)
    end
    float_fields = [fields[i] for (i, t) in enumerate(ftypes) if t <: AbstractFloat]
    kwargs = [:($f = getfield(p, $(QuoteNode(f)))) for f in float_fields]
    quote
        MultiLayer{$FT,N,$P}(; $(kwargs...))
    end
end


has_definedbounds(x) = false
has_definedbounds(x::AbstractLayers) = true

# 多层参数收集：从 S 取字段集和 bounds/units 元数据
@generated function _ml_float_fields(::Type{S}) where {S}
    keep = [f for (f, t) in zip(fieldnames(S), fieldtypes(S)) if t <: AbstractFloat]
    :($(Tuple(keep)))
end

function get_params(x::MultiLayer{FT,N,S}; path=[], with_unit=true) where {FT,N,S}
    fields = _ml_float_fields(S)
    res = map(fields) do field
        value = getproperty(x, field)
        _path = [path..., field]
        bound = bounds(S, field)
        unit = with_unit ? units(S, field) : ""
        [(; path=[_path..., i], name=field,
            value=value[i], type=eltype(value), bound, unit) for i in 1:N]
    end
    filter(r -> !isnothing(r.bound), reduce(vcat, res))
end

# 专属 update!：MultiLayer 的字段（如 :θ_sat）不是真实字段，需走 getproperty
function update!(x::MultiLayer, path::AbstractVector, value::FT; type::Type) where {FT}
    if length(path) >= 2
        field = path[1]
        idx = path[2]
        vec = getproperty(x, field)
        vec[idx] = type(value)
        return
    end
    error("update! on MultiLayer requires path of length >= 2, got $path")
end
