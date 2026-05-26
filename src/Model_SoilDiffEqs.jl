export AbstractRetention, AbstractRetentionLayers
export VanGenuchten, VanGenuchtenLayers, Campbell, CampbellLayers
export HydraulicProfile, _retention_method

export default_hydraulic
export _sync_ksat!
export AbstractThermal, AbstractThermalLayers
export ThermalMain, ThermalMainLayers, ThermalBase, ThermalBaseLayers, ThermalProfile

export SoilModel

abstract type AbstractRetention{T<:Real} end
abstract type AbstractThermal{T<:Real} end

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

# 类型别名
const VanGenuchtenLayers{FT,N} = MultiLayer{FT,N,VanGenuchten{FT}}
const CampbellLayers{FT,N} = MultiLayer{FT,N,Campbell{FT}}

# 抽象 UnionAll 别名（用于 isa/<: 检查和 method dispatch 约束）
const AbstractRetentionLayers{FT,N} = MultiLayer{FT,N,S} where {S<:AbstractRetention{FT}}

_retention_method(::CampbellLayers) = "Campbell"
_retention_method(::VanGenuchtenLayers) = "van_Genuchten"


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

# 从 retention SoA 派生 Ksat 剖面
function KvLayers(retention::MultiLayer{FT,N,S}) where {FT,N,S<:AbstractRetention{FT}}
    KvLayers{FT,N}(; kv=deepcopy(retention.Ksat))
end

_sync_ksat!(kv, layers, dz_cm) = nothing


## Hydraulic Profile
# T 自动从 profile 的 typeof 推断，吸收掉 MultiLayer 的 NT 类型参数
mutable struct HydraulicProfile{FT<:AbstractFloat,N,
    P<:AbstractRetention{FT},
    T<:MultiLayer{FT,N,P},
    K}
    dz_cm::Vector{FT}
    profile::T
    layers::Vector{P}
    kv::K
end

# 外构造器：N 从 profile 类型推断；T/K 从入参 typeof 特化
function HydraulicProfile{FT}(; N::Int=5, dz_cm::Vector{FT}=FT[],
    profile=CampbellLayers{FT,N}(),
    kv=KvLayers{FT,N}()) where {FT<:AbstractFloat}

    Np = length(profile)   # 从 profile 类型参数取 N，忽略 N 入参
    layers = Vector(profile)
    P = eltype(layers)
    HydraulicProfile{FT,Np,P,typeof(profile),typeof(kv)}(dz_cm, profile, layers, kv)
end


## Thermal
@bounds @units @with_kw mutable struct ThermalMain{FT<:AbstractFloat} <: AbstractThermal{FT}
    κ::FT = FT(0.2) | (0.1, 10.0) | "W m-1 K-1"
    cv::FT = FT(2.0e6) | (1.0e6, 5.0e6) | "J m-3 K-1"
end

@bounds @with_kw mutable struct ThermalBase{FT<:AbstractFloat} <: AbstractThermal{FT}
    κ_dry::FT = FT(0.2) | (0.05, 0.5)      # dry soil thermal conductivity [W m-1 K-1]
    ρ_soil::FT = FT(1300.0) | (800.0, 1800.0) # soil bulk density [kg m-3]
    V_SOM::FT = FT(0.02) | (0.0, 0.3)      # organic matter volume fraction [-]
end

const ThermalMainLayers{FT,N} = MultiLayer{FT,N,ThermalMain{FT}}
const ThermalBaseLayers{FT,N} = MultiLayer{FT,N,ThermalBase{FT}}

const AbstractThermalLayers{FT,N} = MultiLayer{FT,N,S} where {S<:AbstractThermal{FT}}


# ThermalProfile可以暂时不使用
mutable struct ThermalProfile{FT<:AbstractFloat,N,
    P<:AbstractThermal{FT},
    T<:MultiLayer{FT,N,P}}
    profile::T
    layers::Vector{P}
end

# 外构造器：N 从 profile 类型推断；T 从入参 typeof 特化
function ThermalProfile{FT}(;
    N::Int=5,
    profile=ThermalMainLayers{FT,N}()) where {FT<:AbstractFloat}

    Np = length(profile)
    layers = Vector(profile)
    P = eltype(layers)
    ThermalProfile{FT,Np,P,typeof(profile)}(profile, layers)
end


##
@with_kw mutable struct SoilModel{FT<:AbstractFloat,H<:HydraulicProfile{FT},T<:ThermalProfile{FT}}
    N::Int = 5                      # number of soil layers
    hydraulic::H = HydraulicProfile{FT}(; N)   # 水力剖面
    thermal::T = ThermalProfile{FT}(; N)       # 热力剖面
end

# 外构造器：不标注 hydraulic/thermal 类型，避免 kwarg 注解 widen 和 T 无法从默认值绑定的问题
function SoilModel{FT}(;
    N::Int=5,
    hydraulic=HydraulicProfile{FT}(; N),
    thermal=ThermalProfile{FT}(; N)) where {FT<:AbstractFloat}

    SoilModel{FT,typeof(hydraulic),typeof(thermal)}(N, hydraulic, thermal)
end
