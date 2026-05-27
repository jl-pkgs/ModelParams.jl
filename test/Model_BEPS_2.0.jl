# Builds on Model_BEPS.jl (ParamVeg, ParamSoilThermal, ParamSoilThermalLayers)
# and Model_SoilDiffEqs.jl (Campbell, CampbellLayers, AbstractRetentionProfile)
include("Model_BEPS_base.jl")

@bounds @with_kw mutable struct ParamBEPS2{FT<:AbstractFloat,
    H<:HydraulicProfile{FT}, T<:ThermalBaseLayers{FT}}

    N::Int = 5
    dz::Vector{FT} = FT[0.05, 0.10, 0.20, 0.40, 1.25]  # 土壤层厚度 [m], BEPS V2023
    r_drainage::FT = Cdouble(0.50) | (0.2, 0.7)  # ? 地表排水速率（地表汇流），可考虑采用曼宁公式

    ψ_min::FT = Cdouble(33.0)  # [m], about 0.10~0.33 MPa开始胁迫点
    alpha::FT = Cdouble(0.4)   # [-], 土壤水限制因子参数，He 2017 JGR-B, Eq. 4

    hydraulic::H = HydraulicProfile{FT,N}()
    thermal::T = ThermalBaseLayers{FT,N}()

    veg::ParamVeg{FT} = ParamVeg{FT}()
end

function ParamBEPS2{FT}(; N::Int=5,
    hydraulic=HydraulicProfile{FT,N}(),
    thermal=ThermalBaseLayers{FT,N}(),
    kwargs...) where {FT<:AbstractFloat}
    ParamBEPS2{FT,typeof(hydraulic),typeof(thermal)}(; N, hydraulic, thermal, kwargs...)
end

##
@testset "Model_BEPS_2.0" begin
    model = ParamBEPS2{Float64}()
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

    # 类型稳定性测试：构造器返回类型须完全 concrete，无自由 _A 等
    @test (@inferred ParamBEPS2{Float64}()) isa ParamBEPS2{Float64}
end

# @code_warntype ParamBEPS2{Float64}(;)
