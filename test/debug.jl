using ModelParams, Parameters
# _sync_ksat!(hydraulic.kv, hydraulic.layers, dz_cm)

N = 4
FT = Float64
p = Campbell(; b=2.0)
profile = Layers(p, N)

kv = KvLayers(profile)
kv = KvExp()
hydraulic = HydraulicProfile{FT,N}(profile, kv)
@time parameters(hydraulic)


##
model = SoilColumn{FT,N}(hydraulic)
@time parameters(model)

##
# @code_warntype kv = KvLayers(profile)
# @code_warntype hydraulic = HydraulicProfile{FT,N}(profile, kv)
# @code_warntype model = SoilColumn{FT,N}(hydraulic)
# @code_warntype parameters(model)
