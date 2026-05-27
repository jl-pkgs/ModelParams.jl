export Retention, Retention_K, Retention_θ, Retention_ψ, Retention_∂θ∂ψ, Retention_∂K∂Se
export Retention_∂K∂θ, Retention_∂ψ∂θ

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
