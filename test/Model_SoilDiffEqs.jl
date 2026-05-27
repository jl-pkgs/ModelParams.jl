using ModelParams, Parameters, Test

_fieldtypes_are_concrete(x) = all(isconcretetype, fieldtypes(typeof(x)))
_struct_is_stable(x) = isconcretetype(typeof(x)) && _fieldtypes_are_concrete(x)

@testset "Model_SoilDiffEqs struct stability" begin
    FT = Float64
    N = 4

    p = Campbell(; b=2.0)
    retention = Layers(p, N)
    kv = @inferred KvLayers(retention)
    hydraulic = HydraulicProfile{FT,N}(retention, kv)
    thermal = ThermalProfile{FT,N}()
    model = SoilModel{FT,N}(hydraulic, thermal)

    @testset "concrete instance layouts" begin
        @test _struct_is_stable(p)
        @test _struct_is_stable(retention)
        @test _struct_is_stable(kv)
        @test _struct_is_stable(hydraulic)
        @test _struct_is_stable(thermal)
        @test _struct_is_stable(model)
    end

    @testset "expected concrete parameterization" begin
        @test retention isa CampbellLayers{FT,N}
        @test kv isa KvLayers{FT,N}
        @test hydraulic isa HydraulicProfile{FT,N,Campbell{FT},typeof(retention),typeof(kv)}
        @test thermal.profile isa ThermalMainLayers{FT,N}
        @test thermal isa ThermalProfile{FT,N,ThermalMain{FT},typeof(thermal.profile)}
        @test model isa SoilModel{FT,N,typeof(hydraulic),typeof(thermal)}
    end

    @testset "nested fields remain concrete" begin
        @test model.hydraulic === hydraulic
        @test model.thermal === thermal
        @test model.hydraulic.profile === retention
        @test model.hydraulic.kv === kv
        @test model.hydraulic.layers isa Vector{Campbell{FT}}
        @test model.thermal.layers isa Vector{ThermalMain{FT}}
    end

    @testset "constructor inference" begin
        @test (@inferred KvLayers(retention)) isa KvLayers{FT,N}
        @test (@inferred HydraulicProfile{FT,N}()) isa HydraulicProfile{FT,N}
        @test (@inferred ThermalProfile{FT,N}()) isa ThermalProfile{FT,N}
        @test (@inferred SoilModel{FT,N}()) isa SoilModel{FT,N}
    end

    params = parameters(model)
    @test :path in propertynames(params)
    @test :bound in propertynames(params)
    @test size(params, 1) > 0
end
