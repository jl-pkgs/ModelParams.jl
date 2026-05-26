export AbstractRetention, AbstractRetentionLayers
export VanGenuchten, VanGenuchtenLayers, Campbell, CampbellLayers
export HydraulicProfile, _retention_method

export default_hydraulic
export _sync_ksat!
export AbstractThermal, AbstractThermalLayers
export ThermalMain, ThermalMainLayers, ThermalBase, ThermalBaseLayers, ThermalProfile

export SoilModel

abstract type AbstractRetention{T<:Real} end
abstract type AbstractRetentionLayers{FT,S} <: AbstractLayers{FT,S} end        # 持水模型层基类

abstract type AbstractThermal{T<:Real} end
abstract type AbstractThermalLayers{FT,S} <: AbstractLayers{FT,S} end          # 热传导模型层基类

@bounds @units @with_kw mutable struct VanGenuchten{T<:Real} <: AbstractRetention{T}
    θ_sat::T = 0.287 | (0.25, 0.50) | "m3 m-3"   # [m3 m-3]
    θ_res::T = 0.075 | (0.03, 0.20) | "m3 m-3"   # [m3 m-3]
    Ksat::T = 34.0 | nothing | "cm h-1"    # [cm h-1]; overridden by kv_profile when set
    α::T = 0.027 | (0.002, 0.300) | "cm-1"       # [cm-1]
    n::T = 3.96 | (1.05, 4.00) | "-"             # [-]
    m::T = 1.0 - 1.0 / n
    ψ_sat::T = 0.0 | nothing | "cm"              # [cm], optional, 
    θ_fc::T = 0.2 | nothing | "m3 m-3"           # field capacity, optional
end

@bounds @units @with_kw mutable struct Campbell{T<:Real} <: AbstractRetention{T}
    θ_sat::T = 0.287 | (0.25, 0.50) | "m3 m-3"
    ψ_sat::T = -10.0 | (-100.0, -5.0) | "cm"
    Ksat::T = 34.0 | nothing | "cm h-1"   # [cm h-1]; overridden by kv_profile when set
    b::T = 4.0 | (2.0, 15.0) | "-"
    θ_fc::T = 0.2 | nothing | "m3 m-3"           # field capacity, optional
end

@make_layers_struct VanGenuchten VanGenuchtenLayers AbstractRetentionLayers
@make_layers_struct Campbell CampbellLayers AbstractRetentionLayers
_retention_method(::VanGenuchtenLayers) = "van_Genuchten"
_retention_method(::CampbellLayers) = "Campbell"


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

function KvLayers(retention::AbstractRetentionLayers{FT}) where {FT<:AbstractFloat}
    KvLayers{FT,length(retention)}(; kv=deepcopy(retention.Ksat))
end

_sync_ksat!(kv, layers, dz_cm) = nothing


@with_kw mutable struct HydraulicProfile{FT<:AbstractFloat,
    S<:AbstractRetentionLayers{FT},
    P<:AbstractRetention{FT},
    K<:Union{AbstractKv,AbstractKvLayers}}
    N::Int = 5
    dz_cm::Vector{FT} = FT[]   # layer thicknesses [cm]; set to enable integral Ksat for exponential profiles

    profile::S = CampbellLayers{FT,N}()       # SoA: 持水模型层参数
    layers::Vector{P} = Vector(profile)       # AoS: 分层实例
    kv::K = KvLayers{FT,N}()                  # Ksat 深度剖面
end

# 外构造器：从 profile/kv 自动推断 S/P/K，免去手写全部类型参数
function HydraulicProfile{FT}(; N::Int=5, dz_cm::Vector{FT}=FT[],
    profile::AbstractRetentionLayers{FT}=CampbellLayers{FT,N}(),
    kv::Union{AbstractKv,AbstractKvLayers}=KvLayers{FT,N}()) where {FT<:AbstractFloat}

    layers = Vector(profile)
    HydraulicProfile{FT,typeof(profile),eltype(layers),typeof(kv)}(;
        N, dz_cm, profile, layers, kv)
end


## Thermal
@bounds @units @with_kw mutable struct ThermalMain{FT<:AbstractFloat} <: AbstractThermal{FT}
    κ::FT = FT(0.2) | (0.1, 10.0) | "W m-1 K-1"
    cv::FT = FT(2.0e6) | (1.0e6, 5.0e6) | "J m-3 K-1"
end
@make_layers_struct ThermalMain ThermalMainLayers AbstractThermalLayers


@bounds @with_kw mutable struct ThermalBase{FT<:AbstractFloat} <: AbstractThermal{FT}
    κ_dry::FT = FT(0.2) | (0.05, 0.5)      # dry soil thermal conductivity [W m-1 K-1]
    ρ_soil::FT = FT(1300.0) | (800.0, 1800.0) # soil bulk density [kg m-3]
    V_SOM::FT = FT(0.02) | (0.0, 0.3)      # organic matter volume fraction [-]
end
@make_layers_struct ThermalBase ThermalBaseLayers AbstractThermalLayers


# ThermalProfile可以暂时不使用
@with_kw mutable struct ThermalProfile{FT<:AbstractFloat,
    S<:AbstractThermalLayers{FT},
    P<:AbstractThermal{FT}}
    N::Int = 5
    profile::S = ThermalMainLayers{FT,N}()
    layers::Vector{P} = Vector(profile)
end

# 外构造器：从 profile 自动推断 S/P
function ThermalProfile{FT}(; N::Int=5,
    profile::AbstractThermalLayers{FT}=ThermalMainLayers{FT,N}()) where {FT<:AbstractFloat}

    layers = Vector(profile)
    ThermalProfile{FT,typeof(profile),eltype(layers)}(; N, profile, layers)
end


##
@with_kw mutable struct SoilModel{FT<:AbstractFloat, H<:HydraulicProfile{FT}, T<:ThermalProfile{FT}}
    N::Int = 5                      # number of soil layers
    hydraulic::H = HydraulicProfile{FT}(; N)   # 水力剖面
    thermal::T = ThermalProfile{FT}(; N)       # 热力剖面
end

# 外构造器：FT 指定时使用默认 hydraulic/thermal
function SoilModel{FT}(; N::Int=5,
    hydraulic::HydraulicProfile{FT}=HydraulicProfile{FT}(; N),
    thermal::ThermalProfile{FT}=ThermalProfile{FT}(; N)) where {FT<:AbstractFloat}

    SoilModel{FT,typeof(hydraulic),typeof(thermal)}(N, hydraulic, thermal)
end
