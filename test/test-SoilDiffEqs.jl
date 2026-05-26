using ModelParams, Parameters, Test

@testset "Model_SoilDiffEqs (src/)" begin
    FT = Float64

    @testset "Retention 类型层次" begin
        @test CampbellLayers{FT,5} <: AbstractRetentionProfile
        @test VanGenuchtenLayers{FT,5} <: AbstractRetentionProfile
        @test Campbell{FT} <: AbstractRetention
        @test VanGenuchten{FT} <: AbstractRetention
    end

    @testset "Vector(AbstractLayers) 展开 SoA → AoS" begin
        layers = CampbellLayers{FT,4}()
        vec = Vector(layers)
        @test vec isa Vector{Campbell{FT}}
        @test length(vec) == 4
        @test vec[1].θ_sat == layers.θ_sat[1]
        @test vec[3].b == layers.b[3]
    end

    @testset "HydraulicProfile: 默认 + 推断 P" begin
        h = HydraulicProfile{FT}()
        @test h.N == 5
        @test h.profile isa CampbellLayers{FT,5}
        @test h.layers isa Vector{Campbell{FT}}
        @test length(h.layers) == 5
        @test h.kv isa KvLayers{FT,5}
    end

    @testset "HydraulicProfile: 替换 retention 模型" begin
        h = HydraulicProfile{FT}(; N=3, profile=VanGenuchtenLayers{FT,3}())
        @test h.profile isa VanGenuchtenLayers{FT,3}
        @test h.layers isa Vector{VanGenuchten{FT}}
        @test length(h.layers) == 3
    end

    @testset "HydraulicProfile: 替换 kv 剖面" begin
        h = HydraulicProfile{FT}(; N=4, kv=KvExpLayers{FT,4}())
        @test h.kv isa KvExpLayers{FT,4}
    end

    @testset "ThermalProfile" begin
        t = ThermalProfile{FT}()
        @test t.N == 5
        @test t.profile isa ThermalMainLayers{FT,5}
        @test t.layers isa Vector{ThermalMain{FT}}
        @test length(t.layers) == 5
        @test t.layers[1].κ == t.profile.κ[1]
    end

    @testset "parameters() 递归收集" begin
        h = HydraulicProfile{FT}()
        df = parameters(h)
        @test :path in propertynames(df)
        @test :bound in propertynames(df)
        # profile 字段下应有 θ_sat、ψ_sat、b（Campbell 有 bound 的字段）
        path_strs = [string.(p) for p in df.path]
        @test any(p -> :θ_sat in p, df.path)
        @test any(p -> :b in p, df.path)
    end

    @testset "update! 修改 profile 中的层参数" begin
        h = HydraulicProfile{FT}(; N=3)
        params = parameters(h)
        path = [:profile, :θ_sat, 2]
        update!(h, [path], [0.42]; params)
        @test h.profile.θ_sat[2] == 0.42
    end

    @testset "Retention method 标识" begin
        @test _retention_method(CampbellLayers{FT,5}()) == "Campbell"
        @test _retention_method(VanGenuchtenLayers{FT,5}()) == "van_Genuchten"
    end
end
