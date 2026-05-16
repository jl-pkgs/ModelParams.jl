module ModelParams


export AbstractModel
export Params, update!
export @bounds, @units, bounds, units
export par_map, @par

using DataFrames

function unlist(list::Vector)
    res = []
    for x in list
        isa(x, Vector) ? append!(res, unlist(x)) : push!(res, x)
    end
    return map(x -> x, res)
end


import FieldMetadata: @metadata

@metadata bounds nothing
@metadata units ""

abstract type AbstractModel{FT} end

include("parallel.jl")
include("ModelParam.jl")
include("GOF.jl")
include("Optim/Optim.jl")


end
