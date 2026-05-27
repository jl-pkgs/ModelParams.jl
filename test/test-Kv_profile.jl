using ModelParams, Test

@testset "KvExp" begin
    kv = KvExp(kv=10.0, f=0.01)

    # At the surface (z=0): Ksat = kv
    @test kv_at_depth(kv, 1, 0.0) ≈ 10.0

    # At depth 100 cm: Ksat = kv * exp(-f*z) = 10 * exp(-1)
    @test kv_at_depth(kv, 3, 100.0) ≈ 10.0 * exp(-1.0)

    # Strictly decreasing with depth
    k1 = kv_at_depth(kv, 1, 10.0)
    k2 = kv_at_depth(kv, 2, 50.0)
    k3 = kv_at_depth(kv, 3, 200.0)
    @test k1 > k2 > k3
end

@testset "KvExpConst" begin
    kv = KvExpConst(kv=10.0, f=0.01, z_exp=50.0)

    # At z < z_exp: same as exponential
    @test kv_at_depth(kv, 1, 20.0) ≈ 10.0 * exp(-0.01 * 20.0)

    # At z_exp: Ksat = kv * exp(-f * z_exp)
    kv_at_zexp = 10.0 * exp(-0.01 * 50.0)
    @test kv_at_depth(kv, 2, 50.0) ≈ kv_at_zexp

    # At z > z_exp: constant (same as at z_exp)
    @test kv_at_depth(kv, 3, 100.0) ≈ kv_at_zexp
    @test kv_at_depth(kv, 4, 300.0) ≈ kv_at_zexp
end

@testset "KvLayers" begin
    kv = KvLayers{Float64,3}(kv=[5.0, 10.0, 2.0])
    @test kv_at_depth(kv, 1, 0.0) ≈ 5.0
    @test kv_at_depth(kv, 2, 10.0) ≈ 10.0
    @test kv_at_depth(kv, 3, 50.0) ≈ 2.0
end

# @testset "KvExpLayers" begin
#     # Upper 2 layers: layered; below: exponential with f[2]=0.01 from kv[2]
#     kv = KvExpLayers{Float64,3}(kv=[5.0, 8.0, 3.0], f=[0.02, 0.01, 0.01])
#     nlayers_kv = 2
#     z_layered_cm = 40.0   # bottom of layer 2

#     # Layer 1 (≤ nlayers_kv): returns kv[1]
#     @test kv_at_depth(kv, 1, 10.0, nlayers_kv, z_layered_cm) ≈ 5.0
#     # Layer 2 (≤ nlayers_kv): returns kv[2]
#     @test kv_at_depth(kv, 2, 30.0, nlayers_kv, z_layered_cm) ≈ 8.0
#     # Layer 3 (> nlayers_kv): exponential from kv[2] with f[2]
#     z3 = 60.0
#     expected = 8.0 * exp(-0.01 * (z3 - z_layered_cm))
#     @test kv_at_depth(kv, 3, z3, nlayers_kv, z_layered_cm) ≈ expected
# end

@testset "kv_layer_ksat — integral formula" begin
    # KvExp: layer-average = kv/f/dz * (exp(-f*z1) - exp(-f*z2))
    kv = KvExp(kv=10.0, f=0.01)
    z1, z2 = 0.0, 20.0
    expected = 10.0 / (0.01 * 20.0) * (exp(-0.01 * z1) - exp(-0.01 * z2))
    @test kv_layer_ksat(kv, 1, z1, z2) ≈ expected

    # Thin layer: degenerates to centre-point
    z1, z2 = 10.0, 10.0 + 1e-10
    @test kv_layer_ksat(kv, 1, z1, z2) ≈ 10.0 * exp(-0.01 * (z1 + z2) / 2)

    # KvExpConst: layer entirely above z_exp
    kvc = KvExpConst(kv=10.0, f=0.01, z_exp=50.0)
    @test kv_layer_ksat(kvc, 1, 0.0, 20.0) ≈
          10.0 / (0.01 * 20.0) * (1.0 - exp(-0.01 * 20.0))

    # Layer entirely below z_exp: constant
    ksat_const = 10.0 * exp(-0.01 * 50.0)
    @test kv_layer_ksat(kvc, 2, 60.0, 80.0) ≈ ksat_const

    # Layer straddling z_exp
    ksat_kvc = kv_layer_ksat(kvc, 1, 40.0, 60.0)
    # First 10 cm: exponential from 40 to 50; last 10 cm: constant
    exp_part = 10.0 / (0.01 * 10.0) * (exp(-0.01 * 40) - exp(-0.01 * 50)) * 10.0 / 20.0
    const_part = ksat_const * 10.0 / 20.0
    @test ksat_kvc ≈ exp_part + const_part

    # _sync_ksat! with dz_cm correctly sets param[i].Ksat via integral
    par = Campbell(θ_sat=0.4, ψ_sat=-10.0, Ksat=10.0, b=4.0)
    kv2 = KvExp(kv=10.0, f=0.01)
    dz = [20.0, 20.0, 20.0]   # [cm]
    # ps = SoilModel(par, 3; kv_profile=kv2, dz_cm=dz)
    # for i in 1:3
    #     z1_cm = (i - 1) * 20.0
    #     z2_cm = i * 20.0
    #     @test ps.param_hydraulic[i].Ksat ≈ kv_layer_ksat(kv2, i, z1_cm, z2_cm)
    # end
end

@testset "KvProfile ModelParams bounds" begin
    kv = KvExp{Float64}()
    # get_opt_info returns (x0, lb, ub, paths) in field-definition order: kv, f
    x0, lb, ub, paths = get_opt_info(kv)
    @test length(x0) == 2
    @test lb[1] ≈ 0.002   # kv lower
    @test ub[1] ≈ 100.0   # kv upper
    @test lb[2] ≈ 0.0     # f lower
    @test ub[2] ≈ 0.1     # f upper

    # KvExpConst has 3 params
    kvc = KvExpConst{Float64}()
    x0c, lbc, ubc, _ = get_opt_info(kvc)
    @test length(x0c) == 3
    @test lbc[3] ≈ 10.0   # z_exp lower
    @test ubc[3] ≈ 500.0  # z_exp upper
end

# @testset "KvProfile layer ModelParams bridge" begin
#     kv = KvExpLayers{Float64,3}(kv=[5.0, 8.0, 3.0], f=[0.02, 0.01, 0.01])
#     @test kv isa AbstractKvLayers

#     params = parameters(kv)
#     @test params.value == [5.0, 8.0, 3.0, 0.02, 0.01, 0.01]
#     @test params.path == [[:kv, 1], [:kv, 2], [:kv, 3], [:f, 1], [:f, 2], [:f, 3]]

#     par = Campbell(θ_sat=0.4, ψ_sat=-10.0, Ksat=10.0, b=4.0)
#     ps = SoilModel(par, 3; kv_profile=kv, nlayers_kv=2, dz_cm=[20.0, 20.0, 20.0])
#     @test ps.kv_profile isa AbstractKvLayers
#     @test get_params(ps, :kv_profile).value == params.value
# end
