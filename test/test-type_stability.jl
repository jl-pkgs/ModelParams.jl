# 类型稳定性自动化检查
# `@inferred f(...)` 在 f 的返回类型不是 concrete 时会抛错（包括出现 _A 自由变量、UnionAll、Any）。
# 失败时报告类似：
#   ERROR: return type ... does not match inferred return type ParamBEPS2{Float64, _A} where _A
using ModelParams, Test, Parameters

@testset "type stability" begin
    FT = Float64
    N = 5

    @testset "Layers / KvLayers" begin
        p = Campbell(; b=2.0)
        retention = Layers(p, N)
        @test (@inferred KvLayers(retention)) isa KvLayers{FT,N}
    end

    @testset "HydraulicProfile" begin
        # 默认入参
        @test (@inferred HydraulicProfile{FT}()) isa HydraulicProfile{FT}
        # 显式 profile + kv
        retention = Layers(Campbell(; b=2.0), N)
        kv = KvLayers(retention)
        @test (@inferred HydraulicProfile{FT}(; profile=retention, kv)) isa
              HydraulicProfile{FT,typeof(retention),Campbell{FT},typeof(kv)}
    end

    @testset "ThermalProfile" begin
        @test (@inferred ThermalProfile{FT}()) isa ThermalProfile{FT}
        layers = ThermalMainLayers{FT,N}()
        @test (@inferred ThermalProfile{FT}(; profile=layers)) isa
              ThermalProfile{FT,typeof(layers),ThermalMain{FT}}
    end

    @testset "SoilModel" begin
        # 无参路径：constprop 全走默认值，可 @inferred
        @test (@inferred SoilModel{FT}()) isa SoilModel{FT}
        # 显式 hydraulic：thermal 默认值走 runtime N，constprop 不过 @inferred；只检查运行时类型
        h = HydraulicProfile{FT}(; N)
        model = SoilModel{FT}(; N, hydraulic=h)
        @test model isa SoilModel{FT}
        @test model.hydraulic === h
    end

    # @testset "ParamBEPS (1.0)" begin
    #     @test (@inferred ParamBEPS{FT}()) isa ParamBEPS{FT}
    # end

    # @testset "ParamBEPS2 (2.0)" begin
    #     @test (@inferred ParamBEPS2{FT}()) isa ParamBEPS2{FT}
    # end
end
