using Random
using Statistics

function geometric_range(x::AbstractMatrix{T}, bound::AbstractVector{T}) where {T<:AbstractFloat}
    # xnstd = std(x, dims=1)[:]                               # standard deviation for each parameter
    gnrng = exp(mean(log.((colMax(x) - colMin(x)) ./ bound))) # normalized geometric range of the parameters
    return gnrng
end

function _callback(nloop, num_evals, fevals)
    bestf = minimum(fevals)
    worstf = maximum(fevals)
    @printf("Iteration = %3d | nEvals = %4d |  Best = %9.5f  |  Worst = %9.5f\n",
        nloop, num_evals, bestf, worstf)
end


"""
SplitMix64 finalizer:
  shifts (30, 27, 31) and multipliers are a validated parameter set,
  used to improve avalanche behavior for nearby input seeds.

避免“相邻种子太像”

# Reference
1. Steele, Lea, Flood (2014), *Fast Splittable Pseudorandom Number Generators* (OOPSLA).
2. Sebastiano Vigna, SplitMix64 / xoshiro reference implementations and notes.
"""
@inline function _splitmix64(x::UInt64)
    z = x + 0x9e3779b97f4a7c15
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end

# 每个 (nloop, igs, iloop) 都有固定且不同的子种子
# Deterministic sub-seed mapping for each (nloop, igs, iloop, tag).
# 0x9e3779b97f4a7c15 is the 64-bit golden-ratio increment commonly used in SplitMix.
@inline function _rng_seed(base_seed::Integer, nloop::Integer, igs::Integer, iloop::Integer=0, tag::Integer=0)
    x = reinterpret(UInt64, Int64(base_seed))
    x ⊻= UInt64(nloop) * 0x9e3779b97f4a7c15
    x ⊻= UInt64(igs) * 0xbf58476d1ce4e5b9
    x ⊻= UInt64(iloop) * 0x94d049bb133111eb
    x ⊻= UInt64(tag) * 0xd2b74407b1ce6e93
    return _splitmix64(x)
end

# 从整数种子创建 RNG, UInt64转成 MersenneTwister实例
# Convert UInt64 sub-seed to a local RNG instance.
@inline _rng_from_seed(seed::UInt64) = Random.MersenneTwister(Int(seed % UInt64(typemax(Int))))


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
