using ModelParams, Parameters, DataFrames, Test

@testset "SoilColumn" begin
    FT = Float64; N = 5
    ps = SoilColumn{FT,N}()

    # parameters() 不受 filter_params 影响（无 dispatch 冲突）
    @test nrow(parameters(ps)) == 30   # 15 Campbell + 5 KvLayers + 10 ThermalMain

    @testset "filter_params mod filter" begin
        @test nrow(filter_params(ps, :hydraulic)) == 20
        @test nrow(filter_params(ps, :thermal))   == 10
        @test nrow(filter_params(ps, :all))        == 30
        @test nrow(filter_params(ps, [:hydraulic, :thermal])) == 30
    end

    @testset "filter_params kwargs" begin
        # list_fix 排除指定字段
        p = filter_params(ps, :hydraulic; list_fix=[:b])
        @test !(:b in p.name) && nrow(p) == 15

        # list_sameLayer 每个名称只保留 1 行
        p = filter_params(ps, :hydraulic; list_sameLayer=[:ψ_sat])
        @test count(==(:ψ_sat), p.name) == 1 && nrow(p) == 16

        # inds 只保留指定层
        p = filter_params(ps, :hydraulic; inds=[1, 2])
        @test nrow(p) == 8
        @test all(x -> x[end] in (1, 2), p.path)
    end

    @testset "update_params!" begin
        # 恒等 round-trip
        ps2 = SoilColumn{FT,N}()
        p = filter_params(ps2, :hydraulic)
        x0 = Float64.(p.value)
        update_params!(ps2, p.path, x0; params=p)
        @test filter_params(ps2, :hydraulic).value ≈ x0

        # list_sameLayer 广播到所有层
        ps3 = SoilColumn{FT,N}()
        ps3.hydraulic.profile.ψ_sat .= FT.(-10.0 .* (1:N))
        p = filter_params(ps3, :hydraulic; list_sameLayer=[:ψ_sat])
        theta = copy(Float64.(p.value))
        theta[findfirst(==(:ψ_sat), p.name)] = FT(-77.0)
        update_params!(ps3, p.path, theta; params=p, list_sameLayer=[:ψ_sat])
        @test all(==(FT(-77.0)), ps3.hydraulic.profile.ψ_sat)
        @test all(l -> l.ψ_sat == FT(-77.0), ps3.hydraulic.layers)
    end

    @testset "ThermalBaseLayers" begin
        ps_tb = SoilColumn{FT,N}(HydraulicProfile{FT,N}(), ThermalProfile{FT,N}(ThermalBaseLayers{FT,N}()))
        @test nrow(parameters(ps_tb)) == 35        # 20 hydraulic + 15 ThermalBase
        p = filter_params(ps_tb, :thermal)
        @test nrow(p) == 15
        @test Set(p.name) == Set([:κ_dry, :ρ_soil, :V_SOM])
    end
end
