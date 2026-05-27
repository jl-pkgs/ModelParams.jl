function cal_se(θ::Real, par::VanGenuchten)
  (; θ_sat, θ_res) = par
  se = (θ - θ_res) / (θ_sat - θ_res)
  return se
end

@testset "Retention_∂K∂θ" begin
  par = VanGenuchten()
  θ1, θ2 = 0.21, 0.22
  se1 = cal_se(θ1, par)
  se2 = cal_se(θ2, par)
  θ = 0.5 * (θ1 + θ2)
  se = 0.5 * (se1 + se2)
  
  K1 = van_Genuchten_K(θ1, par)
  K2 = van_Genuchten_K(θ2, par)
  obs_∂K∂Se = (K2 - K1) / (se2 - se1)
  obs_∂K∂θ = (K2 - K1) / (θ2 - θ1)
  obs_∂K∂Se, obs_∂K∂θ

  sim_∂K∂Se = Retention_∂K∂Se(se, par)
  sim_∂K∂θ = Retention_∂K∂θ(θ, par)

  @test isapprox(sim_∂K∂Se, obs_∂K∂Se, rtol=0.01)
  @test isapprox(sim_∂K∂θ, obs_∂K∂θ, rtol=0.01)
end


@testset "Retention_∂ψ∂θ" begin
  par = VanGenuchten()
  θ1, θ2 = 0.21, 0.22
  ψ1 = Retention_ψ(θ1, par)
  ψ2 = Retention_ψ(θ2, par)  
  obs_∂ψ∂θ = (ψ2 - ψ1) / (θ2 - θ1)

  ψ = 0.5 * (ψ1 + ψ2)
  sim_∂ψ∂θ = Retention_∂ψ∂θ(ψ, par)
  @test isapprox(obs_∂ψ∂θ, sim_∂ψ∂θ, rtol=0.01)
end
