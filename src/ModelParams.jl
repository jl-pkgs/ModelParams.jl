module ModelParams

export AbstractModel
export Params, update!, @bounds, @units, bounds, units
export par_map, @par

using Parameters, DataFrames



function unlist(list::Vector)
    res = []
    for x in list
        isa(x, Vector) ? append!(res, unlist(x)) : push!(res, x)
    end
    return map(x -> x, res)
end


function get_bound(bound::Vector)
    lower = map(x -> x[1], bound)
    upper = map(x -> x[2], bound)
    lower, upper
end

get_bound(params::DataFrame) = get_bound(params.bound)


abstract type AbstractModel{FT} end

include("metadata.jl")

include("parallel.jl")
include("MultiLayer.jl")
include("Params.jl")

include("GOF.jl")
include("Optim/Optim.jl")

include("Kv_Profile.jl")
include("Kv.jl")
include("Model_SoilDiffEqs.jl")
include("SoilColumn.jl")

include("Retention.jl")
include("Retention_Campbell.jl")
include("Retention_van_Genuchten.jl")
include("Retention_ParamTable.jl")
include("Retention_PTF.jl")

include("tridiagonal_solver.jl")

export bounds, units, @bounds, @units



end
