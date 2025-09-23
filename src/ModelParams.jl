module ModelParams


export AbstractModel
export Params, update!
export @bounds, @units

using DataFrames

include("ModelParam.jl")
include("GOF.jl")
include("Optim/Optim.jl")


end
