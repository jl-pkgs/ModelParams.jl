module ModelParams


export AbstractModel
export Params, update!, optimize_params!
export @bounds, @units

include("ModelParam.jl")
include("GOF.jl")
include("Optim/Optim.jl")


end
