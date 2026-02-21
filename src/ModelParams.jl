module ModelParams


export AbstractModel
export Params, update!
export @bounds, @units
export par_map, @par

using DataFrames

include("parallel.jl")
include("ModelParam.jl")
include("GOF.jl")
include("Optim/Optim.jl")


end
