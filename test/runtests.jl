using ModelParams, Test
using Parameters
# 必须要@with_kw，Base.@kwdef报错

include("test-sceua.jl")

@testset "GOF" begin
  @test GOF(1:10, 2:11) ==
        (NSE=0.8787878787878788, R2=1.0, KGE=0.8181818181818181, R=1.0, RMSE=1.0, MAE=1.0, bias=1.0, bias_perc=18.181818181818183, n_valid=10)
end


@testset "ModelParams update!" begin
  include("model_PML.jl")

  FT = Float64
  model = Photosynthesis_Rong2018{FT}()

  params = Params(model)
  # params |> DataFrame

  parnames = [:kQ, :VCmax25, :VPDmin]
  inds = indexin(parnames, params.name)
  paths = params.path[inds]
  
  parvalues = [0.6, 10., 0.8]
  @time update!(model, paths, parvalues; params)

  @test model.kQ == 0.6
  @test model.VCmax25 == 10.
  @test model.watercons.VPDmin == 0.8
end
