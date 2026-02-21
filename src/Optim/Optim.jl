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


ok(msg::String) = printstyled("$msg\n", color=:green)
warn(msg::String) = printstyled("$msg\n", color=:yellow)
bad(msg::String) = printstyled("$msg\n", color=:red)

function show_status(exitflag, maxn)
  if exitflag == ReturnCode.MaxIters
    warn("Optim Stop: maxn=$maxn reached!")
  elseif exitflag == ReturnCode.Success
    ok("Optim Success!")
  elseif exitflag == ReturnCode.Stalled
    ok("Optim Stop: f_reltol reached!")
  elseif exitflag == ReturnCode.Failure
    bad("Optim Failed!")
  end
  # @printf("Search was stopped at trial number: %d \n", num_evals)
  # println("Normalized geometric range = $(num2str(gnrng))")
  # println("The best point has improved in last $(num2str(stagnation_iters)) LOOPS BY $(num2str(criter_change))")
end



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

include("cceua.jl")
include("sceua.jl")

# export SORT, MEAN, MAX, MIN

export sceua, ReturnCode
