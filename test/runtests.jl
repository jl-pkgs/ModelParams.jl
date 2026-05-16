using ModelParams, Test, Parameters
# 必须要@with_kw，Base.@kwdef报错

include("Model_PML.jl")
include("Model_BEPS.jl")
include("Model_SoilDiffEqs.jl")


@testset "GOF" begin
    @test GOF(1:10, 2:11) ==
          (NSE=0.8787878787878788, R2=1.0, KGE=0.8181818181818181, R=1.0, RMSE=1.0, MAE=1.0, bias=1.0, bias_perc=18.181818181818183, n_valid=10)
end

include("test-par_map.jl")
include("test-sceua.jl")
