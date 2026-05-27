include("Model_BEPS_base.jl")

# ParamBEPS2 <: AbstractSoilModel{FT,N}:
#   - thermal 改用 ThermalProfile（统一 SoA + AoS 缓存）
#   - 继承 filter_params / update_params! 无需额外代码
@bounds @with_kw mutable struct ParamBEPS2{FT<:AbstractFloat,N,
    H<:HydraulicProfile{FT,N},T<:ThermalProfile{FT,N}} <: AbstractSoilModel{FT,N}

    dz::Vector{FT} = FT[0.05, 0.10, 0.20, 0.40, 1.25]  # 土壤层厚度 [m]
    r_drainage::FT = 0.50 | (0.2, 0.7) # ? 地表排水速率（地表汇流），可考虑采用曼宁公式

    ψ_min::FT = Cdouble(33.0)  # [m], about 0.10~0.33 MPa开始胁迫点
    alpha::FT = Cdouble(0.4)   # [-], 土壤水限制因子参数，He 2017 JGR-B, Eq. 4

    hydraulic::H
    thermal::T
    veg::ParamVeg{FT} = ParamVeg{FT}()
end

function ParamBEPS2{FT,N}(
    hydraulic=HydraulicProfile{FT,N}(),
    thermal=ThermalProfile{FT,N}(ThermalBaseLayers{FT,N}());
    kwargs...) where {FT<:AbstractFloat,N}
    ParamBEPS2{FT,N,typeof(hydraulic),typeof(thermal)}(; hydraulic, thermal, kwargs...)
end

##
FT = Float64
N = 5
model = ParamBEPS2{FT,N}()
@time parameters(model)
# @code_warntype model = ParamBEPS2{FT,N}()

model = ParamBEPS2{Float64,5}()
params = parameters(model)

opts = [
    (; path=[:r_drainage], name=:r_drainage, value=model.r_drainage),
    (; path=[:veg, :Ω], name=:Ω, value=model.veg.Ω),
    (; path=[:veg, :g1_w], name=:g1_w, value=model.veg.g1_w),
    (; path=[:veg, :g0_w], name=:g0_w, value=model.veg.g0_w),
    (; path=[:veg, :VCmax25], name=:VCmax25, value=model.veg.VCmax25),
] |> DataFrame
paths = opts.path

theta = [0.6, 0.9, 10.0, 0.01, 100.0]
update!(model, paths, theta; params)
par = parameters(model; paths)
@test par.value == theta

##
@testset "Model_BEPS_2.0" begin
    model = ParamBEPS2{FT,N}()

    # filter_params 共用同一实现（来自 AbstractSoilModel / SoilColumn.jl）
    @test nrow(filter_params(model, :hydraulic)) == 20   # 15 Campbell + 5 KvLayers
    @test nrow(filter_params(model, :thermal)) == 15   # ThermalBase: κ_dry, ρ_soil, V_SOM × 5
    @test all(x -> x[1] === :veg, filter_params(model, :veg).path)

    # list_fix / list_sameLayer 同样适用
    p = filter_params(model, :hydraulic; list_fix=[:b], list_sameLayer=[:ψ_sat])
    @test !(:b in p.name)
    @test count(==(:ψ_sat), p.name) == 1

    # update_params! 共用同一实现（AbstractSoilModel dispatch）
    p = filter_params(model, :hydraulic)
    x0 = Float64.(p.value)
    update_params!(model, p.path, x0; params=p)
    @test filter_params(model, :hydraulic).value ≈ x0

    # AoS 缓存同步
    @test model.hydraulic.layers[1].θ_sat == model.hydraulic.profile.θ_sat[1]
    @test model.thermal.layers[1].κ_dry == model.thermal.profile.κ_dry[1]

    # 类型稳定性
    @test (@inferred ParamBEPS2{FT,N}()) isa ParamBEPS2{FT,N}
end
