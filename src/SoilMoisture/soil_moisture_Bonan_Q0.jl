function soil_moisture_Q0!(
  θ::V, ψ::V, θ_prev::V, ψ_prev::V, ∂θ∂ψ::V, K::V, K₊ₕ::V, tri::TriSolver{FT},
  ps::HydraulicProfile, sink::V, Q0::FT;
  ibeg::Int, N::Int, Δz_cm::V, Δz₊ₕ_cm::V, dt::FT) where {FT<:Real,V<:AbstractVector{FT}}

  (; u, a, b, c, d, e, f) = tri
  ψ_next = u
  dt = dt / 3600  # [s] -> [h]

  θ_prev[ibeg:N] .= θ[ibeg:N]
  ψ_prev[ibeg:N] .= ψ[ibeg:N]

  cal_ψ!(ψ, ps, θ; N, ibeg)
  cal_K!(K, K₊ₕ, ps, θ; N, ibeg, Δz=Δz_cm)
  cal_∂θ∂ψ!(∂θ∂ψ, ps, ψ; N, ibeg)

  @inbounds for i = ibeg:N
    if i == ibeg
      a[i] = 0
      c[i] = -K₊ₕ[i] / Δz₊ₕ_cm[i]
      b[i] = ∂θ∂ψ[i] * Δz_cm[i] / (0.5 * dt) - c[i]
      d[i] = ∂θ∂ψ[i] * Δz_cm[i] / (0.5 * dt) * ψ[i] - K₊ₕ[i] - Q0
    elseif i < N
      a[i] = -K₊ₕ[i-1] / Δz₊ₕ_cm[i-1]
      c[i] = -K₊ₕ[i] / Δz₊ₕ_cm[i]
      b[i] = ∂θ∂ψ[i] * Δz_cm[i] / (0.5 * dt) - a[i] - c[i]
      d[i] = ∂θ∂ψ[i] * Δz_cm[i] / (0.5 * dt) * ψ[i] + K₊ₕ[i-1] - K₊ₕ[i]
    elseif i == N
      a[i] = -K₊ₕ[N-1] / Δz₊ₕ_cm[N-1]
      c[i] = 0
      b[i] = ∂θ∂ψ[i] * Δz_cm[i] / (0.5 * dt) - a[i] - c[i]
      d[i] = ∂θ∂ψ[i] * Δz_cm[i] / (0.5 * dt) * ψ[i] + K₊ₕ[N-1] - K[i]
    end
    d[i] -= sink[i]
  end

  tridiagonal_solver!(a, b, c, d, e, f, ψ_next; ibeg, N)

  cal_θ!(θ, ps, ψ_next; N, ibeg)
  cal_K!(K, K₊ₕ, ps, θ; N, ibeg, Δz=Δz_cm)
  cal_∂θ∂ψ!(∂θ∂ψ, ps, ψ_next; N, ibeg)

  @inbounds for i = ibeg:N
    if i == ibeg
      a[i] = 0
      c[i] = -K₊ₕ[i] / (2 * Δz₊ₕ_cm[i])
      b[i] = ∂θ∂ψ[i] * Δz_cm[i] / dt - c[i]
      d[i] = ∂θ∂ψ[i] * Δz_cm[i] / dt * ψ[i] - Q0 +
             c[i] * (ψ[i] - ψ[i+1]) - K₊ₕ[i]
    elseif i < N
      a[i] = -K₊ₕ[i-1] / (2 * Δz₊ₕ_cm[i-1])
      c[i] = -K₊ₕ[i] / (2 * Δz₊ₕ_cm[i])
      b[i] = ∂θ∂ψ[i] * Δz_cm[i] / dt - a[i] - c[i]
      d[i] = ∂θ∂ψ[i] * Δz_cm[i] / dt * ψ[i] - a[i] * (ψ[i-1] - ψ[i]) +
             c[i] * (ψ[i] - ψ[i+1]) + K₊ₕ[i-1] - K₊ₕ[i]
    else
      i == N
      a[i] = -K₊ₕ[i-1] / (2 * Δz₊ₕ_cm[i-1])
      c[i] = 0
      b[i] = ∂θ∂ψ[i] * Δz_cm[i] / dt - a[i] - c[i]
      d[i] = ∂θ∂ψ[i] * Δz_cm[i] / dt * ψ[i] - a[i] * (ψ[i-1] - ψ[i]) + K₊ₕ[i-1] - K[i]
    end
    d[i] -= sink[i]
  end

  tridiagonal_solver!(a, b, c, d, e, f, ψ; ibeg, N)
  cal_θ!(θ, ps, ψ; N, ibeg)

  QN = -K[N]
  dθ = 0
  for i = ibeg:N
    dθ += (θ[i] - θ_prev[i]) * Δz_cm[i]
  end

  err = dθ - (QN - Q0) * dt
  Q0, QN, dθ, err
end
