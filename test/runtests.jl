using ModelParams, Test
using Parameters
# 必须要@with_kw，Base.@kwdef报错

include("test-sceua.jl")


@testset "GOF" begin
  @test GOF(1:10, 2:11) ==
        (NSE=0.8787878787878788, R2=1.0, KGE=0.8181818181818181, R=1.0, RMSE=1.0, MAE=1.0, bias=1.0, bias_perc=18.181818181818183, n_valid=10)
end


# 水分限制-GPP
abstract type AbstractWaterConsGPPModel{FT} <: AbstractModel{FT} end

# 光合
abstract type AbstractPhotosynthesisModel{FT} <: AbstractModel{FT} end


# # 气孔导度
# abstract type AbstractStomatalModel{FT} <: AbstractModel{FT} end
@bounds @units @with_kw mutable struct β_GPP_Zhang2019{FT} <: AbstractWaterConsGPPModel{FT}
  ## water constraint
  "parameter to constrain `gc`"
  VPDmin::FT = 0.9 | (0.65, 1.5) | "kPa"

  "parameter to constrain `gc`"
  VPDmax::FT = 4.0 | (3.50, 6.5) | "kPa"
end


@bounds @units @with_kw mutable struct Photosynthesis_Rong2018{FT} <: AbstractPhotosynthesisModel{FT}
  "initial slope of the light response curve to assimilation rate, (i.e., quantum efficiency)"
  α::FT = 0.06 | (0.01, 0.10) | "μmol CO2 [μmol PAR]⁻¹"

  "initial slope of the CO2 response curve to assimilation rate, (i.e., carboxylation efficiency)"
  η::FT = 0.04 | (0.01, 0.07) | "μmol m⁻² s⁻¹ [μmol m⁻² s⁻¹]⁻¹"

  "carbon saturated rate of photosynthesis at 25 °C"
  VCmax25::FT = 50.00 | (5.00, 120.00) | "μmol m⁻² s⁻¹"

  "photoperiod constraint"
  d_pc::FT = 2.0 | (0.0, 5.0) | "-"

  "extinction coefficients for visible radiation" # 植被光合参数
  kQ::FT = 0.45 | (0.10, 1.0) | "-"

  watercons::AbstractWaterConsGPPModel{FT} = β_GPP_Zhang2019{FT}()
end


@testset "ModelParams update!" begin
  FT = Float64
  model = Photosynthesis_Rong2018{FT}()

  params = Params(model)
  # params |> DataFrame

  parnames = [:kQ, :VCmax25, :VPDmin]
  parvalues = [0.6, 10., 0.8]
  @time update!(model, parnames, parvalues; params)

  @test model.kQ == 0.6
  @test model.VCmax25 == 10.
  @test model.watercons.VPDmin == 0.8
end


@testset "Parameter Optimization" begin
  FT = Float64
  model = Photosynthesis_Rong2018{FT}()
  
  # Test optimization function exists and can be called
  @test hasmethod(optimize_params!, (AbstractModel, Function, Vector{Symbol}))
  
  # Simple objective function to minimize squared distance from targets
  objective_fn = function(m)
    target_α = 0.08
    target_kQ = 0.5
    return (m.α - target_α)^2 + (m.kQ - target_kQ)^2
  end
  
  initial_cost = objective_fn(model)
  
  # Test optimization with manual bounds
  param_names = [:α, :kQ]
  manual_bounds = [(0.01, 0.10), (0.10, 1.0)]
  
  optimal_params, optimal_cost, exitflag = optimize_params!(
    model, objective_fn, param_names; 
    bounds=manual_bounds, maxn=200, verbose=false
  )
  
  # Basic checks
  @test length(optimal_params) == 2
  @test optimal_cost < initial_cost
  @test optimal_params[1].name == :α
  @test optimal_params[2].name == :kQ
  
  # Test error handling
  @test_throws ErrorException optimize_params!(model, objective_fn, [:nonexistent])
end
