using ModelParams, Test, Parameters
# 必须要@with_kw，Base.@kwdef报错

include("test-Kv_profile.jl")

include("test-SoilDiffEqs.jl")
include("test-SoilColumn.jl")
include("Model_SoilDiffEqs.jl")

##
include("Model_PML.jl")
include("Model_BEPS_base.jl")
include("Model_BEPS_1.0.jl")
include("Model_BEPS_2.0.jl")

include("test-type_stability.jl")

##
@testset "GOF" begin
    @test GOF(1:10, 2:11) ==
          (NSE=0.8787878787878788, R2=1.0, KGE=0.8181818181818181, R=1.0, RMSE=1.0, MAE=1.0, bias=1.0, bias_perc=18.181818181818183, n_valid=10)
end

include("test-par_map.jl")
include("test-sceua.jl")
