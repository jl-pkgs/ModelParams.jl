using ModelParams, Parameters, DataFrames, Test
import ModelParams: bounds, units, @bounds, @units

abstract type AbstractSoilParam{T<:Real} end


@bounds @units @with_kw mutable struct VanGenuchten{T<:Real} <: AbstractSoilParam{T}
    θ_sat::T = 0.287 | (0.25, 0.50) | "m3 m-3"   # [m3 m-3]
    θ_res::T = 0.075 | (0.03, 0.20) | "m3 m-3"     # [m3 m-3]
    Ksat::T = 34.0 | (0.002, 60.0) | "cm h-1"        # [cm h-1]
    α::T = 0.027 | (0.002, 0.300) | "cm-1"        # [cm-1]
    n::T = 3.96 | (1.05, 4.00) | "-"         # [-]
    m::T = 1.0 - 1.0 / n
end

@bounds @units @with_kw mutable struct Campbell{T<:Real} <: AbstractSoilParam{T}
    θ_sat::T = 0.287 | (0.25, 0.50) | "m3 m-3"
    ψ_sat::T = -10.0 | (-100.0, -5.0) | "cm"
    Ksat::T = 34.0 | (0.002, 100.0) | "cm h-1"
    b::T = 4.0 | (3.0, 15.0) | "-"
end

@make_layers_struct VanGenuchten
@make_layers_struct Campbell

const ParamSoilHydraulic{FT,N} = Union{VanGenuchtenLayers{FT,N},CampbellLayers{FT,N}} where {FT<:AbstractFloat,N}

@with_kw mutable struct ParamSoil{FT<:AbstractFloat,N}
    hydraulic::ParamSoilHydraulic{FT,N} = VanGenuchtenLayers{FT,N}()
    # thermal::ParamSoilThermal{FT} = ParamSoilThermal{FT}()
end


@generated function Layers(p::P, N::Int) where {FT, P<:AbstractSoilParam{FT}}
    LayersName = Symbol(P.name.name, :Layers)
    kwargs = [:($(f) = fill(p.$f, N)) for f in fieldnames(P)]
    quote
        $LayersName{$FT,N}($(kwargs...))
    end
end

function ParamSoil(p::P, N::Int) where {FT, P<:AbstractSoilParam{FT}}
    ParamSoil{FT,N}(hydraulic=Layers(p, N))
end


##
model = ParamSoil{Float64,1}()
model = ParamSoil{Float64,4}()
parameters(model)

##
p = Campbell(; b=2.0)
model = ParamSoil(p, 4)
parameters(model)
