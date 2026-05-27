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
        @test (@inferred HydraulicProfile{FT,N}()) isa HydraulicProfile{FT,N}
        # 显式 profile + kv
        retention = Layers(Campbell(; b=2.0), N)
        kv = KvLayers(retention)
        @test (@inferred HydraulicProfile{FT,N}(retention, kv)) isa
              HydraulicProfile{FT,N,Campbell{FT},typeof(retention),typeof(kv)}
    end

    @testset "ThermalProfile" begin
        @test (@inferred ThermalProfile{FT,N}()) isa ThermalProfile{FT,N}
        layers = ThermalMainLayers{FT,N}()
        @test (@inferred ThermalProfile{FT,N}(layers)) isa
              ThermalProfile{FT,N,ThermalMain{FT},typeof(layers)}
    end

    @testset "SoilModel" begin
        # 无参路径：constprop 全走默认值，可 @inferred
        @test (@inferred SoilModel{FT,N}()) isa SoilModel{FT,N}
        # 显式 hydraulic：N 作为类型参数传入，thermal 从 hydraulic 的类型参数 N 派生
        h = HydraulicProfile{FT,N}()
        model = @inferred SoilModel{FT,N}(h)
        @test model isa SoilModel{FT,N}
        @test model.hydraulic === h
        @test model.thermal isa ThermalProfile{FT,N}
    end

    # @testset "ParamBEPS (1.0)" begin
    #     @test (@inferred ParamBEPS{FT}()) isa ParamBEPS{FT}
    # end

    # @testset "ParamBEPS2 (2.0)" begin
    #     @test (@inferred ParamBEPS2{FT}()) isa ParamBEPS2{FT}
    # end
end
