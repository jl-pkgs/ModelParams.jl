export build_params, ParamThermal
using ModelParams, Parameters
import ModelParams: bounds, units, @bounds, @units, @make_layers_struct


##
@bounds @with_kw mutable struct SoilModel{FT<:AbstractFloat}
    N::Int = 5   # number of soil layers
    dz_cm::Vector{FT} = FT[]                                # layer thicknesses [cm]; set to enable integral Ksat for exponential profiles

    hydraulic::HydraulicProfile{FT,<:AbstractRetention{FT}} =
        HydraulicProfile{FT}(; N, profile=Layers(default_hydraulic(FT), N))
    thermal::ThermalProfile{FT} =
        ThermalProfile{FT}(; N, profile=Layers(ParamThermal{FT}(), N))
end

function SoilModel(retention::AbstractRetentionLayers{FT}, N::Int;
    thermal=Layers(ParamThermal{FT}(), N), 
    kv=default_kv_profile(retention),
    dz_cm::AbstractVector=FT[]) where {FT<:AbstractFloat}

    hydraulic_profile = HydraulicProfile{FT}(; N, profile=retention, kv)
    thermal_profile = ThermalProfile{FT}(; N, profile=thermal)

    _sync_ksat!(hydraulic_profile.kv, hydraulic_profile.layers, dz_cm)
    SoilModel{FT}(N, dz_cm, hydraulic_profile, thermal_profile)
end

# 默认开启的是多层参数
function SoilModel(p::AbstractRetention{FT}, N::Int; kw...) where {FT<:AbstractFloat}
    retention = Layers(p, N)
    SoilModel(retention, N; kw...)
end

##
# p = Campbell(; b=2.0)
# model = SoilModel(p, 4)
# parameters(model)
