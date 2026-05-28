"""
    Campbell(ψ::T, par::Campbell{T})

Campbell (1974) relationships

# Arguments
+ `ψ`: Matric potential, cm
+ `param`
  - `ψ_sat`: Matric potential at saturation, [cm]
  - `θ_sat`: Volumetric water content at saturation
  - `b`    : Exponent
  - `Ksat`: Hydraulic conductivity at saturation [cm h-1]

# TODO: 核对变量的单位

# Examples
```julia
θ_sat = 0.25
ψ_sat = -25.0
b = 0.2
Ksat = 3.4e-03
θ, K, ∂θ∂ψ = Campbell(ψ, ψ_sat, θ_sat, Ksat, b)
```
"""
@inline function Campbell(ψ::T, par::Campbell{T}) where {T<:Real}
  (; ψ_sat, θ_sat, Ksat, b) = par
  if ψ <= ψ_sat
    ratio = ψ / ψ_sat
    Se = ratio^(-one(T) / b)
    θ = θ_sat * Se
    K = Ksat * Se^(2b + 3)
    ∂θ∂ψ = -θ_sat / (b * ψ_sat) * Se / ratio
    return θ, K, ∂θ∂ψ
  end
  θ_sat, Ksat, zero(T)
end

"""
    Campbell_θ(ψ, ψ_sat, θ_sat, b)
"""
@inline @fastmath function Campbell_θ(ψ::T, par::Campbell{T}) where {T<:Real}
  (; ψ_sat, θ_sat, b) = par
  ψ <= ψ_sat ? θ_sat * (ψ / ψ_sat)^(-1 / b) : θ_sat
end

"""
    Campbell_K(θ, θ_sat, Ksat, b)
"""
@inline @fastmath function Campbell_K(θ::T, par::Campbell{T}) where {T<:Real}
  (; θ_sat, Ksat, b) = par
  Se = clamp(θ / θ_sat, T(0.01), T(1.0))
  Ksat * Se^(2b + 3)
end


@inline @fastmath function Campbell_∂θ∂ψ(ψ::T, par::Campbell{T}) where {T<:Real}
  (; ψ_sat, θ_sat, b) = par
  ψ <= ψ_sat ? -θ_sat / (b * ψ_sat) * (ψ / ψ_sat)^(-1 / b - 1) : T(0.0)
end

@inline Campbell_∂ψ∂θ(ψ::T, par::Campbell{T}) where {T<:Real} = T(1.0) / Campbell_∂θ∂ψ(ψ, par)


# Campbell 1974, Bonan 2019 Table 8.2
@inline @fastmath function Campbell_ψ(θ::T, par::Campbell{T}; ψ_min=T(-1e7)) where {T<:Real}
  (; θ_sat, ψ_sat, b) = par
  Se = clamp(θ / θ_sat, T(0.01), T(1.0))
  ψ = ψ_sat * Se^(-b)
  return max(ψ, ψ_min) # ψ为负值
end

@inline @fastmath function Campbell_ψ_Se(Se::T, par::Campbell{T}; ψ_min=T(-1e7)) where {T<:Real}
  (; ψ_sat, b) = par
  Se = clamp(Se, T(0.01), T(1.0))
  ψ = ψ_sat * Se^(-b)
  return max(ψ, ψ_min) # ψ为负值
end

@inline @fastmath function Campbell_∂K∂Se(Se::T, par::Campbell{T}) where {T<:Real}
  (; Ksat, b) = par
  Ksat * (2b + 3) * (Se^(2b + 2))
end

@inline @fastmath function Campbell_∂K∂θ(θ::T, par::Campbell{T}) where {T<:Real}
  (; Ksat, b) = par
  Se = clamp(θ / par.θ_sat, T(0.01), T(1.0))
  Campbell_∂K∂Se(Se, par) / (par.θ_sat)
end

export Campbell, Campbell_ψ, Campbell_θ, Campbell_K,
  Campbell_∂θ∂ψ, Campbell_∂ψ∂θ, Campbell_∂K∂Se, Campbell_ψ_Se
