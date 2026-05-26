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
