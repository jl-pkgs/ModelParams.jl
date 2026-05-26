export AbstractKvProfile, AbstractKv, AbstractKvLayers
export Kv, KvLayers, KvExp, KvExpLayers, KvExpConst, KvExpPiecewise

abstract type AbstractKvProfile{T<:Real} end
abstract type AbstractKv{FT<:Real} <: AbstractKvProfile{FT} end


"""Per-layer Ksat (scalar stub). Use KvLayers for multi-layer instances."""
@bounds @units @with_kw mutable struct Kv{T<:Real} <: AbstractKv{T}
    kv::T = 34.0 | (0.002, 60.0) | "cm h-1"
end


"""Exponential decline: Ksat(z) = kv · exp(−f · z)"""
@bounds @units @with_kw mutable struct KvExp{T<:Real} <: AbstractKv{T}
    kv::T = 34.0 | (0.002, 100.0) | "cm h-1"   # surface Ksat
    f::T = 0.01 | (0.0, 0.1) | "cm-1"     # depth-decay coefficient
end


"""Exponential + constant: exponential for z < z_exp, then constant below."""
@bounds @units @with_kw mutable struct KvExpConst{T<:Real} <: AbstractKv{T}
    kv::T = 34.0 | (0.002, 100.0) | "cm h-1"
    f::T = 0.01 | (0.0, 0.1) | "cm-1"
    z_exp::T = 100.0 | (10.0, 500.0) | "cm"    # depth at which exponential stops
end


# 类型别名：用 MultiLayer 统一表达
const KvLayers{FT,N}        = MultiLayer{FT,N,Kv{FT}}
const KvExpLayers{FT,N}     = MultiLayer{FT,N,KvExp{FT}}
const KvExpPiecewise{FT,N}  = MultiLayer{FT,N,KvExpConst{FT}}

# 抽象 UnionAll 别名：所有 Kv 系列的 MultiLayer
const AbstractKvLayers{FT,N} = MultiLayer{FT,N,S,NT} where {S<:AbstractKv{FT},NT<:NamedTuple}
