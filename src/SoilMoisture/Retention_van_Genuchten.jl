# 2.5x faster power method
"Faster method for exponentiation"
pow(x, y) = x^y
# @fastmath pow(x::Real, y::Real) = exp(y * log(x))

"""
    van_Genuchten(ψ::T, par::VanGenuchten{T})

van Genuchten (1980) relationships

# Arguments
+ `ψ`: Matric potential
+ `param`
  - `θ_res`       : Residual water content
  - `θ_sat`       : Volumetric water content at saturation
  - `Ksat`        : Hydraulic conductivity at saturation [cm h-1]
  - `α`           : Inverse of the air entry potential (cm-1)
  - `n`           : Pore-size distribution index
  - `m`           : Exponent
  - `soil_texture`: Soil texture flag

# Examples
```julia
# Haverkamp et al. (1977): sand
param = (soil_texture = 1, 
  θ_res = 0.075, θ_sat = 0.287, 
  α = 0.027, n = 3.96, m = 1, Ksat = 34)

# Haverkamp et al. (1977): Yolo light clay
param = (soil_texture=2, 
  θ_res = 0.124, θ_sat = 0.495,
  α = 0.026, n = 1.43, m = 1 - 1 / 1.43,
  Ksat = 0.0443)
```
"""
@inline function van_Genuchten(ψ::T, par::VanGenuchten{T}) where {T<:Real}
  (; θ_res, θ_sat, Ksat, α, n, m) = par

  if ψ <= T(-1e7)
    return θ_res, zero(T), van_Genuchten_∂θ∂ψ(ψ, par)
  elseif ψ > zero(T)
    return θ_sat, Ksat, zero(T)
  end

  x = -α * ψ
  xnm1 = x^(n - one(T))
  xn = x * xnm1
  den = one(T) + xn
  Se = den^(-m)
  θ = θ_res + (θ_sat - θ_res) * Se

  diff = xn / den
  K = Ksat * sqrt(Se) * (one(T) - diff^m)^2
  ∂θ∂ψ = α * m * n * (θ_sat - θ_res) * xnm1 * Se / den
  θ, K, ∂θ∂ψ
end

# @fastmath 
function van_Genuchten_θ(ψ::T, par::VanGenuchten{T}; ψ_min::T=T(-1e7)) where {T<:Real}
  (; θ_res, θ_sat, α, n, m) = par
  ψ <= ψ_min && return θ_res

  # Effective saturation (Se) for specified matric potential (ψ)
  Se = ψ <= 0 ? (1 + (α * abs(ψ))^n)^-m : 1.0
  Se = clamp(Se, T(0.0), T(1.0))

  # Volumetric soil moisture (θ) for specified matric potential (ψ)
  return θ_res + (θ_sat - θ_res) * Se # θ
end

# @fastmath 
function van_Genuchten_K(θ::T, par::VanGenuchten{T}) where {T<:Real}
  (; θ_res, θ_sat, Ksat, m) = par
  Se::T = (θ - θ_res) / (θ_sat - θ_res)
  Se = clamp(Se, T(0.0), T(1.0))

  diff::T = (1 - Se^(1 / m))
  K::T = Se < 1 ? Ksat * sqrt(Se) * (1 - diff^m)^2 : Ksat
  return K
end

# @fastmath 
@inline function van_Genuchten_∂θ∂ψ(ψ::T, par::VanGenuchten{T})::T where {T<:Real}
  (; θ_res, θ_sat, α, n, m) = par
  if ψ <= zero(T)
    # num = α * m * n * (θ_sat - θ_res) * (α * abs(ψ))^(n - 1)
    # den = (1 + (α * abs(ψ))^n)^(m + 1)
    x = -α * ψ
    term1 = x^(n - one(T))
    num = α * m * n * (θ_sat - θ_res) * term1
    den = (one(T) + term1 * x)^(m + one(T))
    ∂θ∂ψ = num / den
  else
    ∂θ∂ψ = zero(T)
  end
  return ∂θ∂ψ
end

van_Genuchten_∂ψ∂θ(ψ::T, par::VanGenuchten{T}) where {T<:Real} = T(1.0) / van_Genuchten_∂θ∂ψ(ψ, par)


# @fastmath 
# ψmin = -1e7cm, CLM5, Eq. 7.53
function van_Genuchten_ψ(θ::T, par::VanGenuchten{T}; ψ_min=T(-1e7)) where {T<:Real}
  (; θ_res, θ_sat, α, n, m) = par
  if θ <= θ_res
    return ψ_min # Return a very high negative number indicating very dry conditions
  elseif θ >= θ_sat
    return T(0.0)   # Saturated condition, psi is zero
  else
    Se = (θ - θ_res) / (θ_sat - θ_res)
    Se = clamp(Se, T(0.0), T(1.0))
    ψ = -1 / α * pow(pow(1.0 / Se, (1 / m)) - 1, 1 / n)
    return max(ψ, ψ_min) # Ensure the returned value does not go below ψ_min
  end
end

function van_Genuchten_ψ_Se(Se::T, par::VanGenuchten{T}; ψ_min=T(-1e7)) where {T<:Real}
  (; α, n, m) = par
  Se <= 0.0 && return ψ_min
  Se >= 1.0 && return T(0.0)
  ψ = -1 / α * pow(pow(1.0 / Se, (1 / m)) - 1, 1 / n)
  return max(ψ, ψ_min) # Ensure the returned value does not go below ψ_min
end


@inline function van_Genuchten_∂K∂Se(Se::T, par::VanGenuchten{T}) where {T<:Real}
  (; Ksat, m) = par
  f = 1 - (1 - Se^(1 / m))^m
  term1 = f^2 / (2 * sqrt(Se))
  term2 = 2 * Se^(1 / m - 1 / 2) * f / ((1 - Se^(1 / m))^(1 - m))
  return Ksat * (term1 + term2)
end

@inline function van_Genuchten_∂K∂θ(θ::T, par::VanGenuchten{T}) where {T<:Real}
  (; θ_res, θ_sat) = par
  Se = (θ - θ_res) / (θ_sat - θ_res)
  return van_Genuchten_∂K∂Se(Se, par) / (θ_sat - θ_res)
end


export van_Genuchten, van_Genuchten_θ, van_Genuchten_K, van_Genuchten_ψ,
  van_Genuchten_∂θ∂ψ, van_Genuchten_∂ψ∂θ, van_Genuchten_∂K∂Se,
  van_Genuchten_ψ_Se, 
  van_Genuchten_∂K∂θ

# Special case for:
# - `soil_type = 1`: Haverkamp et al. (1977) sand
# - `soil_type = 2`: Yolo light clay

# if soil_type == 1
#   K = Ksat * 1.175e6 / (1.175e6 + abs(ψ)^4.74)
# elseif soil_type == 2
#   K = Ksat * 124.6 / (124.6 + abs(ψ)^1.77)
# end
