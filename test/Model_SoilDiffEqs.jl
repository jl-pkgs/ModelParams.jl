export build_params, ParamThermal
export AbstractKvProfile, AbstractKv, AbstractKvLayers
export Kv, KvLayers, KvExp, KvExpLayers, KvExpConst, KvExpPiecewise

using ModelParams, Parameters
import ModelParams: bounds, units, @bounds, @units, @make_layers_struct

abstract type AbstractSoilParam{T<:Real} end
abstract type AbstractKvProfile{T<:Real} end
abstract type AbstractKv{FT<:Real} <: AbstractKvProfile{FT} end
abstract type AbstractKvLayers{FT,S} <: AbstractLayers{FT,S} end

function Base.getindex(x::AbstractLayers{FT,S}, i::Int) where {FT,S}
    kw = (; (name => getfield(x, name)[i] for name in fieldnames(typeof(x)))...)
    return S{FT}(; kw...)
end

Base.length(x::AbstractLayers) = length(getfield(x, first(fieldnames(typeof(x)))))

function build_params(x::AbstractLayers{FT,S}, N::Int) where {FT,S}
    Np = length(x)
    Type = typeof(x)

    map(i -> begin
            Np == 1 && (i = 1) # 如果所有层参数相同，则始终取第一个
            kw = (name => getfield(x, name)[i] for name in fieldnames(Type)) |> NamedTuple
            S{FT}(; kw...)
        end, 1:N)
end

##
@bounds @units @with_kw mutable struct VanGenuchten{T<:Real} <: AbstractSoilParam{T}
    θ_sat::T = 0.287 | (0.25, 0.50) | "m3 m-3"   # [m3 m-3]
    θ_res::T = 0.075 | (0.03, 0.20) | "m3 m-3"   # [m3 m-3]
    Ksat::T = 34.0 | nothing | "cm h-1"    # [cm h-1]; overridden by kv_profile when set
    α::T = 0.027 | (0.002, 0.300) | "cm-1"       # [cm-1]
    n::T = 3.96 | (1.05, 4.00) | "-"             # [-]
    m::T = 1.0 - 1.0 / n
    ψ_sat::T = 0.0 | nothing | "cm"              # [cm], optional, 
    θ_fc::T = 0.2 | nothing | "m3 m-3"           # field capacity, optional
end

@bounds @units @with_kw mutable struct Campbell{T<:Real} <: AbstractSoilParam{T}
    θ_sat::T = 0.287 | (0.25, 0.50) | "m3 m-3"
    ψ_sat::T = -10.0 | (-100.0, -5.0) | "cm"
    Ksat::T = 34.0 | nothing | "cm h-1"   # [cm h-1]; overridden by kv_profile when set
    b::T = 4.0 | (2.0, 15.0) | "-"
    θ_fc::T = 0.2 | nothing | "m3 m-3"           # field capacity, optional
end

@make_layers_struct VanGenuchten
@make_layers_struct Campbell
_retention_method(::VanGenuchtenLayers) = "van_Genuchten"
_retention_method(::CampbellLayers) = "Campbell"

const RetentionLayers{FT,N} = Union{VanGenuchtenLayers{FT,N},CampbellLayers{FT,N}}

function default_hydraulic(::Type{FT}, method_retention::String="van_Genuchten") where {FT<:AbstractFloat}
    if method_retention == "van_Genuchten"
        p = VanGenuchten{FT}(; θ_sat=0.4, θ_res=0.1, Ksat=2.0, α=0.01, n=2.0)
    elseif method_retention == "Campbell"
        p = Campbell{FT}(; θ_sat=0.4, ψ_sat=-10.0, Ksat=2.0, b=4.0)
    else
        error("Unknown method_retention: $method_retention")
    end
    return p
end


##
@bounds @units @with_kw mutable struct ParamThermal{FT<:AbstractFloat}
    κ::FT = FT(0.2) | (0.1, 10.0) | "W m-1 K-1"
    cv::FT = FT(2.0e6) | (1.0e6, 5.0e6) | "J m-3 K-1"
end
@make_layers_struct ParamThermal ParamThermalLayers

## ThermalProfile: 绑定 properties（SoA）与 layers（AoS）
@with_kw mutable struct ThermalProfile{FT<:AbstractFloat}
    properties::ParamThermalLayers{FT}        # SoA: 热参数层 (κ, cv)
    layers::Vector{ParamThermal{FT}}          # AoS: 分层实例
end

function ThermalProfile(properties::ParamThermalLayers{FT}, N::Int) where {FT<:AbstractFloat}
    Np = length(properties)
    if Np == 1
        properties = Layers(properties[1], N)
    elseif Np != N
        error("thermal layers length ($Np) must be 1 or match soil layers ($N).")
    end
    layers = build_params(properties, N)
    ThermalProfile{FT}(properties, layers)
end

## Kv profiles

"""Per-layer Ksat (scalar stub). Use KvLayers for multi-layer instances."""
@bounds @units @with_kw mutable struct Kv{T<:Real} <: AbstractKv{T}
    kv::T = 34.0 | (0.002, 60.0) | "cm h-1"
end
@make_layers_struct Kv KvLayers AbstractKvLayers


"""Exponential decline: Ksat(z) = kv · exp(−f · z)"""
@bounds @units @with_kw mutable struct KvExp{T<:Real} <: AbstractKv{T}
    kv::T = 34.0 | (0.002, 100.0) | "cm h-1"   # surface Ksat
    f::T = 0.01 | (0.0, 0.1) | "cm-1"     # depth-decay coefficient
end
@make_layers_struct KvExp KvExpLayers AbstractKvLayers


"""Exponential + constant: exponential for z < z_exp, then constant below."""
@bounds @units @with_kw mutable struct KvExpConst{T<:Real} <: AbstractKv{T}
    kv::T = 34.0 | (0.002, 100.0) | "cm h-1"
    f::T = 0.01 | (0.0, 0.1) | "cm-1"
    z_exp::T = 100.0 | (10.0, 500.0) | "cm"    # depth at which exponential stops
end
@make_layers_struct KvExpConst KvExpPiecewise AbstractKvLayers

function default_kv_profile(retention::RetentionLayers{FT}, N::Int) where {FT<:AbstractFloat}
    layers = build_params(retention, N)
    KvLayers{FT,N}(; kv=FT[p.Ksat for p in layers])
end

_sync_ksat!(kv, layers, dz_cm) = nothing

## HydraulicProfile: 绑定 retention（SoA）、layers（AoS）、kv（Ksat 深度剖面）
@with_kw mutable struct HydraulicProfile{FT<:AbstractFloat,P<:AbstractSoilParam{FT}}
    retention::RetentionLayers{FT}                       # SoA: 持水模型层参数
    layers::Vector{P}                                     # AoS: 分层实例
    kv::Union{AbstractKv,AbstractKvLayers}                # Ksat 深度剖面
end

function HydraulicProfile(retention::RetentionLayers{FT}, N::Int; kv=nothing) where {FT<:AbstractFloat}
    layers = build_params(retention, N)
    isnothing(kv) && (kv = default_kv_profile(retention, N))
    HydraulicProfile{FT,eltype(layers)}(retention, layers, kv)
end

##
@bounds @with_kw mutable struct SoilModel{FT<:AbstractFloat,P<:AbstractSoilParam{FT}}
    N::Int = 5   # number of soil layers
    Np::Int = N  # number of parameter layers
    dz_cm::Vector{FT} = FT[]                                # layer thicknesses [cm]; set to enable integral Ksat for exponential profiles

    hydraulic::HydraulicProfile{FT,P} = HydraulicProfile(Layers(default_hydraulic(FT), Np), N)
    thermal::ThermalProfile{FT} = ThermalProfile(Layers(ParamThermal{FT}(), Np), N)
end

function SoilModel(retention::RetentionLayers{FT}, N::Int;
    thermal=nothing, kv=nothing,
    dz_cm::AbstractVector=FT[]) where {FT<:AbstractFloat}
    Np = length(retention)
    isnothing(thermal) && (thermal = Layers(ParamThermal{FT}(), Np))

    dz_cm_vec = FT.(dz_cm)
    hydraulic_profile = HydraulicProfile(retention, N; kv)
    thermal_profile = ThermalProfile(thermal, N)
    P = eltype(hydraulic_profile.layers)

    _sync_ksat!(hydraulic_profile.kv, hydraulic_profile.layers, dz_cm_vec)
    SoilModel{FT,P}(N, Np, dz_cm_vec, hydraulic_profile, thermal_profile)
end

# 默认开启的是多层参数
function SoilModel(p::AbstractSoilParam{FT}, N::Int, Np::Int=N; kw...) where {FT<:AbstractFloat}
    retention = Layers(p, Np)
    SoilModel(retention, N; kw...)
end

##
# p = Campbell(; b=2.0)
# model = SoilModel(p, 4)
# parameters(model)
