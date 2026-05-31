export Retention, Retention_K, Retention_θ, Retention_ψ, Retention_∂θ∂ψ, Retention_∂K∂Se
export Retention_∂K∂θ, Retention_∂ψ∂θ
export mean_arithmetic, mean_harmonic
export cal_θ!, cal_∂θ∂ψ!, cal_ψ!, cal_K!, cal_θKCap!, cal_K_CLM5!
export ψ0_param_index, Init_ψ0

# 多重派发(runtime-dispatch)可能会导致速度变慢
Retention(ψ::T, par::Campbell{T}) where {T<:Real} = Campbell(ψ, par)
Retention_K(θ::T, par::Campbell{T}) where {T<:Real} = Campbell_K(θ, par)
Retention_θ(ψ::T, par::Campbell{T}) where {T<:Real} = Campbell_θ(ψ, par)
Retention_ψ(θ::T, par::Campbell{T}) where {T<:Real} = Campbell_ψ(θ, par)
Retention_ψ_Se(Se::T, par::Campbell{T}) where {T<:Real} = Campbell_ψ_Se(Se, par)
Retention_∂K∂θ(θ::T, par::Campbell{T}) where {T<:Real} = Campbell_∂K∂θ(θ, par)
Retention_∂θ∂ψ(ψ::T, par::Campbell{T}) where {T<:Real} = Campbell_∂θ∂ψ(ψ, par)
Retention_∂ψ∂θ(ψ::T, par::Campbell{T}) where {T<:Real} = Campbell_∂ψ∂θ(ψ, par)
Retention_∂K∂Se(Se::T, par::Campbell{T}) where {T<:Real} = Campbell_∂K∂Se(Se, par)


Retention(ψ::T, par::VanGenuchten{T}) where {T<:Real} = van_Genuchten(ψ, par)
Retention_θ(ψ::T, par::VanGenuchten{T}) where {T<:Real} = van_Genuchten_θ(ψ, par)
Retention_K(θ::T, par::VanGenuchten{T}) where {T<:Real} = van_Genuchten_K(θ, par)
Retention_ψ(θ::T, par::VanGenuchten{T}) where {T<:Real} = van_Genuchten_ψ(θ, par)
Retention_ψ_Se(Se::T, par::VanGenuchten{T}) where {T<:Real} = van_Genuchten_ψ_Se(Se, par)
Retention_∂θ∂ψ(ψ::T, par::VanGenuchten{T}) where {T<:Real} = van_Genuchten_∂θ∂ψ(ψ, par)
Retention_∂ψ∂θ(ψ::T, par::VanGenuchten{T}) where {T<:Real} = van_Genuchten_∂ψ∂θ(ψ, par)
Retention_∂K∂Se(Se::T, par::VanGenuchten{T}) where {T<:Real} = van_Genuchten_∂K∂Se(Se, par)
Retention_∂K∂θ(θ::T, par::VanGenuchten{T}) where {T<:Real} = van_Genuchten_∂K∂θ(θ, par)


Retention(ψ::T; par::AbstractRetention{T}) where {T<:Real} = Retention(ψ, par)
Retention_K(θ::T; par::AbstractRetention{T}) where {T<:Real} = Retention_K(θ, par)
Retention_θ(ψ::T; par::AbstractRetention{T}) where {T<:Real} = Retention_θ(ψ, par)
Retention_ψ(θ::T; par::AbstractRetention{T}) where {T<:Real} = Retention_ψ(θ, par)
Retention_ψ_Se(Se::T; par::AbstractRetention{T}) where {T<:Real} = Retention_ψ_Se(Se, par)
Retention_∂θ∂ψ(ψ::T; par::AbstractRetention{T}) where {T<:Real} = Retention_∂θ∂ψ(ψ, par)
Retention_∂ψ∂θ(ψ::T; par::AbstractRetention{T}) where {T<:Real} = Retention_∂ψ∂θ(ψ, par)
Retention_∂K∂Se(Se::T; par::AbstractRetention{T}) where {T<:Real} = Retention_∂K∂Se(Se, par)

mean_arithmetic(K1::T, K2::T, d1::T, d2::T) where {T<:Real} = (K1 * d1 + K2 * d2) / (d1 + d2)
mean_harmonic(K1::T, K2::T, d1::T, d2::T) where {T<:Real} = K1 * K2 * (d1 + d2) / (K1 * d2 + K2 * d1)

# ─── Primary API: output first, ps second, N/ibeg/Δz as kwargs ────────────────
function cal_θ!(θ::V, ps::HydraulicProfile, ψ::V;
  N::Int, ibeg::Int=1) where {T<:AbstractFloat,V<:AbstractVector{T}}

  param = ps.layers
  @inbounds for i = ibeg:N
    θ[i] = Retention_θ(ψ[i], param[i])
  end
end

function cal_∂θ∂ψ!(∂θ∂ψ::V, ps::HydraulicProfile, ψ::V;
  N::Int, ibeg::Int=1) where {T<:AbstractFloat,V<:AbstractVector{T}}

  param = ps.layers
  @inbounds for i = ibeg:N
    ∂θ∂ψ[i] = Retention_∂θ∂ψ(ψ[i], param[i])
  end
end

function cal_ψ!(ψ::V, ps::HydraulicProfile, θ::V;
  N::Int, ibeg::Int=1) where {T<:AbstractFloat,V<:AbstractVector{T}}

  param = ps.layers
  i0 = max(ibeg - 1, 1)
  @inbounds for i = i0:N
    ψ[i] = Retention_ψ(θ[i], param[i])
  end
end

# Requires ps.hydraulic.dz_cm to be set so that _sync_ksat! has pre-baked
# depth-integrated K_sat into ps.layers[i].K_sat. Never call with an empty dz_cm
# HydraulicProfile when kv is not KvLayers.
function cal_K!(K::V, K₊ₕ::V, ps::HydraulicProfile, θ::V;
  N::Int, ibeg::Int=1, Δz::V) where {T<:AbstractFloat,V<:AbstractVector{T}}

  param = ps.layers
  i0 = max(ibeg - 1, 1)
  @inbounds for i = i0:N
    K[i] = Retention_K(θ[i], param[i])
  end
  @inbounds for i = i0:N-1
    K₊ₕ[i] = mean_arithmetic(K[i], K[i+1], Δz[i], Δz[i+1])
  end
  K₊ₕ[N] = K[N]
end

function cal_θKCap!(θ::V, K::V, K₊ₕ::V, ∂θ∂ψ::V,
  ps::HydraulicProfile, ψ::V;
  N::Int, ibeg::Int=1, Δz::V) where {T<:AbstractFloat,V<:AbstractVector{T}}
  param = ps.layers
  i0 = max(ibeg - 1, 1)
  i0 < ibeg && (K[i0] = Retention_K(θ[i0], param[i0]))
  @inbounds for i = ibeg:N
    θ[i], K[i], ∂θ∂ψ[i] = Retention(ψ[i], param[i])
  end
  @inbounds for i = i0:N-1
    K₊ₕ[i] = mean_arithmetic(K[i], K[i+1], Δz[i], Δz[i+1])
  end
  K₊ₕ[N] = K[N]
  return nothing
end

function cal_K_CLM5!(K::V, K₊ₕ::V,
  ps::HydraulicProfile, θ::V;
  N::Int, ibeg::Int=1, Δz_cm::V) where {T<:AbstractFloat,V<:AbstractVector{T}}

  param = ps.layers
  i0 = max(ibeg - 1, 1)
  @inbounds for i = i0:N
    K[i] = Retention_K(θ[i], param[i])
  end

  K[N+1] = Retention_K(θ[N], param[N])
  @inbounds for i = i0:N
    K₊ₕ[i] = mean_arithmetic(K[i], K[i+1], Δz_cm[i], Δz_cm[i+1])
  end
  return K₊ₕ
end

ψ0_param_index(ibeg::Int, ψ0_location::Symbol) =
  ψ0_location === :center ? max(ibeg - 1, 1) : ibeg

function Init_ψ0(ps::HydraulicProfile, θ::T;
  ibeg::Int, ψ0_location::Symbol) where {T<:Real}

  return Retention_ψ(θ, ps.layers[ψ0_param_index(ibeg, ψ0_location)])
end
