using ModelParams, Test, Parameters
# 必须要@with_kw，Base.@kwdef报错

include("Model_SoilDiffEqs.jl")

FT = Float64
N = 5
p = VanGenuchten{FT}(; θ_sat=0.4, θ_res=0.1, Ksat=2.0, α=0.01, n=2.0)
model = SoilModel(p, 1)
model = SoilModel(p, N)
# model = SoilModel{Float64}(; N=1)
# model = SoilModel{Float64}(; N=5)
parameters(model)

##
include("Model_PML.jl")
include("Model_BEPS_base.jl")
include("Model_BEPS_1.0.jl")
include("Model_BEPS_2.0.jl")

##
@testset "GOF" begin
    @test GOF(1:10, 2:11) ==
          (NSE=0.8787878787878788, R2=1.0, KGE=0.8181818181818181, R=1.0, RMSE=1.0, MAE=1.0, bias=1.0, bias_perc=18.181818181818183, n_valid=10)
end

include("test-par_map.jl")
include("test-sceua.jl")
