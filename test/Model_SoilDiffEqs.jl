using ModelParams, Parameters
# _sync_ksat!(hydraulic.kv, hydraulic.layers, dz_cm)

N = 4
FT = Float64
p = Campbell(; b=2.0)
retention = Layers(p, N)

kv = KvLayers(retention)
hydraulic = HydraulicProfile{FT}(; profile=retention, kv)
model = SoilModel{FT}(; N, hydraulic)
parameters(model)


# @code_warntype kv = KvLayers(retention)
# @code_warntype hydraulic = HydraulicProfile{FT}(; profile=retention, kv)
# @code_warntype model = SoilModel{FT}(; N, hydraulic)
# @code_warntype parameters(model)
