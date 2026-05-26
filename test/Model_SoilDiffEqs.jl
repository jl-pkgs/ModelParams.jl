using ModelParams, Parameters
import ModelParams: bounds, units, @bounds, @units, @make_layers_struct


@bounds @with_kw mutable struct SoilModel{FT<:AbstractFloat}
    N::Int = 5   # number of soil layers
    hydraulic::HydraulicProfile{FT,<:AbstractRetention{FT}} = HydraulicProfile{FT}(; N)
    thermal::ThermalProfile{FT} = ThermalProfile{FT}(; N)
end
# _sync_ksat!(hydraulic.kv, hydraulic.layers, dz_cm)

##
N = 4
FT = Float64
p = Campbell(; b=2.0)
retention = Layers(p, N)

kv = KvLayers(retention)
hydraulic = HydraulicProfile{FT}(; profile=retention, kv)

model = SoilModel(; N, hydraulic)
parameters(model)
