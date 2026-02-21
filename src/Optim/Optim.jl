# export SORT, MEAN, MAX, MIN
export sceua, ReturnCode


import Statistics: std, mean;
import Printf: @sprintf, @printf
using Random


disp = println;
num2str = string;

module ReturnCode
@enum T begin
    Default
    Success
    MaxIters
    Stalled
    Failure
end
end

# 兜底机制
sanitize(v::FT) where {FT<:Real} = isfinite(v) ? v : FT(Inf)

function SORT(x::Vector{<:Real})
    idx = sortperm(x)
    x[idx], idx
end

function colMax(x, dims=1)
    maximum(x, dims=dims)[:]
end

function colMin(x, dims=1)
    minimum(x, dims=dims)[:]
end

# function MEAN(x, dims)
#   mean(x, dims=dims)[:]
# end
# MEAN(x) = MEAN(x, 1)

## RANDOM ----------------------------------------------------------------------
set_seed(seed=1) = Random.seed!(seed)
mrand() = rand()
mrand(n) = rand(n)

# using MATLAB
# set_seed(seed=1) = mat"rand('seed', $seed);"
# mrand() = mat"rand();"
# mrand(n) = mat"rand($n, 1);"
# include("MATLAB_helper.jl")

include("helper.jl")
include("cceua.jl")
include("sceua.jl")
