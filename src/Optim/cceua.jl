"""
    cceua(fn, s, sf, bl, bu, icall)

Generate a new point in a simplex

# Example
```julia
snew, fnew, icall = cceua(fn, s, sf, bl, bu, icall)
```
"""
function cceua(fn, s::AbstractMatrix{FT}, sf::AbstractVector{FT},
    bl::Vector{FT}, bu::Vector{FT}, icall, args...;
    rng=Random.default_rng(), kw...) where {FT<:Real}

    nps, n_param = size(s)
    n = nps
    # m = nopt
    alpha = FT(1.0)
    beta = FT(0.5)

    # Assign the best & worst points:
    # sb = s[1, :]
    # fb = sf[1]
    sw = s[n, :]
    fw = sf[n]

    # Compute the centroid of the simplex excluding the worst point:
    ce = mean(s[1:n-1, :], dims=1)[:]
    snew = ce .+ alpha * (ce .- sw)  # Attempt a reflection point

    # Check if is outside the bounds:
    ibound = 0

    s1 = snew - bl
    idx = findall(s1 .< 0)
    !isempty(idx) && (ibound = 1)

    s1 = bu - snew
    idx = findall(s1 .< 0)
    !isempty(idx) && (ibound = 2)

    if ibound >= 1
        snew = bl + rand(rng, FT, n_param) .* (bu - bl)
    end

    fnew = fn(snew, args...; kw...)
    icall += 1

    # Reflection failed; now attempt a contraction point:
    if fnew .> fw
        snew = sw .+ beta * (ce .- sw)
        fnew = fn(snew, args...; kw...)
        icall += 1

        # Both reflection & contraction have failed; attempt a random point
        if fnew .> fw
            snew = bl + rand(rng, FT, n_param) .* (bu - bl)
            fnew = fn(snew, args...; kw...)
            icall += 1
        end
    end
    snew, fnew, icall
end
