function soil_moisture!(
  θ::V, ψ::V, θ_prev::V, ψ_prev::V, ∂θ∂ψ::V, K::V, K₊ₕ::V, tri::TriSolver{FT},
  ps::HydraulicProfile, sink::V, ψ0::FT;
  ibeg::Int, N::Int, Δz_cm::V, Δz₊ₕ_cm::V, dt::FT,
  ψ0_boundary::NamedTuple{(:i0, :dz0₊ₕ),Tuple{Int,FT}}, ψ0_location=:boundary,
  debug::Bool=true) where {FT<:Real,V<:AbstractVector{FT}}

  (; u, a, b, c, d, e, f) = tri
  ψ_next = u
  dt = dt / 3600 # [s] -> [h]

  for i in 1:N  # backup
    θ_prev[i] = θ[i]
    ψ_prev[i] = ψ[i]
  end

  # cal_ψ!(ψ, ps, θ; N, ibeg)
  # cal_K!(K, K₊ₕ, ps, θ; N, ibeg, Δz=Δz_cm)
  # cal_∂θ∂ψ!(∂θ∂ψ, ps, ψ; N, ibeg)
  cal_θKCap!(θ, K, K₊ₕ, ∂θ∂ψ, ps, ψ; N, ibeg, Δz=Δz_cm)

  (; i0, dz0₊ₕ) = ψ0_boundary # TODO: improve this interface
  K0₊ₕ = (ψ0_location === :boundary || ibeg == 1) ? K[ibeg] : K₊ₕ[i0]

  dt_half = 0.5 * dt
  # first round:
  @inbounds for i = ibeg:N
    if i == ibeg
      a[i] = 0.0
      c[i] = -K₊ₕ[i] / Δz₊ₕ_cm[i]
      b[i] = ∂θ∂ψ[i] * Δz_cm[i] / dt_half + K0₊ₕ / dz0₊ₕ - c[i]
      d[i] = ∂θ∂ψ[i] * Δz_cm[i] / dt_half * ψ[i] + K0₊ₕ / dz0₊ₕ * ψ0 + K0₊ₕ - K₊ₕ[i]
    elseif i < N
      a[i] = -K₊ₕ[i-1] / Δz₊ₕ_cm[i-1]
      c[i] = -K₊ₕ[i] / Δz₊ₕ_cm[i]
      b[i] = ∂θ∂ψ[i] * Δz_cm[i] / dt_half - a[i] - c[i]
      d[i] = ∂θ∂ψ[i] * Δz_cm[i] / dt_half * ψ[i] + K₊ₕ[i-1] - K₊ₕ[i]
    elseif i == N
      a[i] = -K₊ₕ[N-1] / Δz₊ₕ_cm[N-1]
      c[i] = 0.0
      b[i] = ∂θ∂ψ[i] * Δz_cm[i] / dt_half - a[i] - c[i]
      d[i] = ∂θ∂ψ[i] * Δz_cm[i] / dt_half * ψ[i] + K₊ₕ[N-1] - K[i]
    end
    d[i] -= sink[i]
  end
  tridiagonal_solver!(a, b, c, d, e, f, ψ_next; ibeg, N)

  ## update: θ, K and ∂θ∂ψ
  # cal_θ!(θ, ps, ψ_next; N, ibeg)
  # cal_K!(K, K₊ₕ, ps, θ; N, ibeg, Δz_cm)
  # cal_∂θ∂ψ!(∂θ∂ψ, ps, ψ_next; N, ibeg)
  cal_θKCap!(θ, K, K₊ₕ, ∂θ∂ψ, ps, ψ_next; N, ibeg, Δz=Δz_cm)
  K0₊ₕ = (ψ0_location === :boundary || ibeg == 1) ? K[ibeg] : K₊ₕ[i0]

  ## second round: in half step
  @inbounds for i = ibeg:N
    if i == ibeg
      a[i] = 0
      c[i] = -K₊ₕ[i] / (2 * Δz₊ₕ_cm[i])
      b[i] = ∂θ∂ψ[i] * Δz_cm[i] / dt - c[i] + K0₊ₕ / (2 * dz0₊ₕ)
      d[i] = ∂θ∂ψ[i] * Δz_cm[i] / dt * ψ[i] +
             K0₊ₕ / (2 * dz0₊ₕ) * (2ψ0 - ψ[i]) +
             c[i] * (ψ[i] - ψ[i+1]) + K0₊ₕ - K₊ₕ[i]
    elseif i < N
      a[i] = -K₊ₕ[i-1] / (2 * Δz₊ₕ_cm[i-1])
      c[i] = -K₊ₕ[i] / (2 * Δz₊ₕ_cm[i])
      b[i] = ∂θ∂ψ[i] * Δz_cm[i] / dt - a[i] - c[i]
      d[i] = ∂θ∂ψ[i] * Δz_cm[i] / dt * ψ[i] - a[i] * (ψ[i-1] - ψ[i]) +
             c[i] * (ψ[i] - ψ[i+1]) + K₊ₕ[i-1] - K₊ₕ[i]
    elseif i == N
      a[i] = -K₊ₕ[i-1] / (2 * Δz₊ₕ_cm[i-1])
      c[i] = 0
      b[i] = ∂θ∂ψ[i] * Δz_cm[i] / dt - a[i] - c[i]
      d[i] = ∂θ∂ψ[i] * Δz_cm[i] / dt * ψ[i] - a[i] * (ψ[i-1] - ψ[i]) + K₊ₕ[i-1] - K[i]
    end
    d[i] -= sink[i]
  end
  tridiagonal_solver!(a, b, c, d, e, f, ψ; ibeg, N)
  # cal_θ!(soil, ψ)

  ## Check water balance
  Q0 = -K0₊ₕ / (2 * dz0₊ₕ) * ((ψ0 - ψ_prev[ibeg]) + (ψ0 - ψ[ibeg])) - K0₊ₕ # 两个时刻的
  QN = -K[N]

  dθ = 0
  for i = 1:N
    dθ += (θ[i] - θ_prev[i]) * Δz_cm[i]
  end

  err = dθ - (QN - Q0) * dt
  debug && return Q0, QN, dθ, err
  return nothing
end
