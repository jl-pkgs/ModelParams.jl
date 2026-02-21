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


function sceua(fn::Function, x0::Vector{FT}, bl::Vector{FT}, bu::Vector{FT}, args...;
  verbose=false, parallel=true,
  maxn=1000, kstop=5, f_reltol=0.0001, x_reltol=0.0001, n_complex=5,
  seed=1, include_initial=1, kw...) where {FT<:Real}

  exitflag = ReturnCode.Default
  set_seed(seed)

  # Initialize SCE parameters:
  n_param = length(x0)
  npg = 2 * n_param + 1
  nps = n_param + 1
  nspl = npg
  n_popu = npg * n_complex

  bound = bu - bl
  # Create an initial population to fill array x[npt,nopt]:

  x = zeros(FT, n_popu, n_param)
  for i = 1:n_popu
    x[i, :] = bl + rand(n_param) .* bound
  end
  include_initial == 1 && (x[1, :] = x0)

  _fn(i, _args...; _kw...) = fn(x[i, :], _args...; _kw...) |> sanitize
  xf = FT.(par_map(_fn, 1:n_popu, args...; kw..., parallel=true, use_deepcopy=true))
  # xf = zeros(FT, npt)
  # for i = 1:npt
  #   xf[i] = fn(x[i, :]) # nopt
  # end

  # 判定标准
  num_evals = n_popu
  gnrng = geometric_range(x, bound)
  criter_change = 1e+5

  xf, idx = SORT(xf)
  x = x[idx, :]

  nloop = 0
  bestx = x[1, :]
  bestf = xf[1]
  @printf("Iteration = %3d, nEvals = %3d, Best Cost = %.5f\n", nloop, num_evals, bestf)

  # Begin evolution loops:
  lpos = 1
  criter = []

  while num_evals .< maxn && gnrng > x_reltol && criter_change .> f_reltol
    nloop = nloop + 1
    # Loop on complexes [sub-populations]
    for igs = 1:n_complex
      # Partition the population into complexes [sub-populations]
      k1 = 1:npg
      k2 = (k1 .- 1) .* n_complex .+ igs
      cx = x[k2, :]
      cf = xf[k2]

      # Evolve sub-population igs for nspl steps:
      lcs = zeros(Int, nps)
      for loop = 1:nspl
        # Select simplex by sampling the complex according to a linear
        # probability distribution
        lcs[1] = 1
        for k3 = 2:nps
          for iter = 1:1000
            lpos = 1 + floor(npg + 0.5 - sqrt((npg + 0.5)^2 - npg * (npg + 1) * rand()))
            idx = findall(lcs[1:k3-1] .== lpos)
            isempty(idx) && break
          end
          lcs[k3] = lpos
        end
        lcs = sort(lcs)
        # Construct the simplex:
        s = zeros(FT, nps, n_param)
        s = cx[lcs, :]
        sf = cf[lcs]

        snew, fnew, num_evals = cceua(fn, s, sf, bl, bu, num_evals, args...; kw...)

        # Replace the simplex into the complex()
        cx[lcs[end], :] = snew
        cf[lcs[end]] = fnew
        # Sort the complex()
        cf, idx = SORT(cf)
        cx = cx[idx, :]
        # End of Inner Loop for Competitive Evolution of Simplexes
      end
      # Replace the complex back into the population
      x[k2, :] = cx[k1, :]
      xf[k2] = cf[k1]
      # End of Loop on Complex Evolution
    end
    # Shuffled the complexes
    xf, idx = SORT(xf)
    x = x[idx, :]

    # Record the best & worst points
    bestx = x[1, :]
    bestf = xf[1]

    gnrng = geometric_range(x, bound)

    @printf("Iteration = %3d, nEvals = %3d, Best Cost = %.5f\n", nloop, num_evals, bestf)

    # Check for convergency
    push!(criter, bestf)
    # criter = [criter bestf]'
    if (nloop >= kstop)
      criter_change = abs(criter[nloop] - criter[nloop-kstop+1])
      criter_change = criter_change / mean(abs.(criter[nloop-kstop+1:nloop]))
    end
    # End of the Outer Loops
  end

  if exitflag == ReturnCode.Default
    if num_evals >= maxn
      exitflag = ReturnCode.MaxIters
    elseif gnrng <= x_reltol
      exitflag = ReturnCode.Success
    elseif criter_change <= f_reltol
      exitflag = ReturnCode.Stalled
    else
      exitflag = ReturnCode.Failure
    end
  end
  verbose && show_status(exitflag, maxn)

  bestx, bestf, exitflag
end
