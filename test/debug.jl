using ModelParams, Parameters
# _sync_ksat!(hydraulic.kv, hydraulic.layers, dz_cm)

N = 4
FT = Float64
p = Campbell(; b=2.0)
profile = Layers(p, N)

kv = KvLayers(profile)
hydraulic = HydraulicProfile{FT,N}(profile, kv)

# model = SoilModel{FT,N}(; hydraulic)
# parameters(model)

@code_warntype kv = KvLayers(profile)
@code_warntype hydraulic = HydraulicProfile{FT,N}(profile, kv)
@code_warntype model = SoilModel{FT,N}(hydraulic)
# @code_warntype parameters(model)
