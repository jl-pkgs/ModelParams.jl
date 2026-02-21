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
  verbose=false,
  maxn=1000, kstop=5, f_reltol=0.0001, x_reltol=0.0001, n_complex=5, 
  seed=1, iniflg=1, kw...) where {FT<:Real}

  exitflag = -1
  set_seed(seed)

  # Initialize SCE parameters:
  nopt = length(x0)
  npg = 2 * nopt + 1
  nps = nopt + 1
  nspl = npg
  npt = npg * n_complex

  bound = bu - bl
  # Create an initial population to fill array x[npt,nopt]:

  x = zeros(FT, npt, nopt)
  for i = 1:npt
    # x[i, :] = bl + rand(nopt) .* bound
    x[i, :] = bl + mrand(nopt) .* bound
  end

  iniflg == 1 && (x[1, :] = x0)

  icall = 0
  xf = zeros(FT, npt)
  for i = 1:npt
    xf[i] = fn(x[i, :]) # nopt
    icall += 1
  end

  # Sort the population in order of increasing function values
  xf, idx = SORT(xf)
  x = x[idx, :]

  # Record the best & worst points
  nloop = 0
  bestx = x[1, :]
  bestf = xf[1]
  @printf("Iteration = %3d, nEvals = %3d, Best Cost = %.5f\n", nloop, icall, bestf)

  gnrng = geometric_range(x, bound)

  # disp("The Initial Loop: 0")
  # disp(["BESTF  : ' num2str(bestf), ' ' 'BESTX  : [' num2str(bestx), ']"])
  # disp(["WORSTF : ' num2str(worstf), ' ' 'WORSTX : [' num2str(worstx), ']"])
  # nloop = 0

  # Check for convergency
  if icall >= maxn
    disp("*** OPTIMIZATION SEARCH TERMINATED BECAUSE THE LIMIT")
    disp("ON THE MAXIMUM NUMBER OF TRIALS $maxn HAS BEEN EXCEEDED.")
    println("SEARCH WAS STOPPED AT TRIAL NUMBER: $icall OF THE INITIAL LOOP!")
  end

  if gnrng .< x_reltol
    exitflag = 1
    disp("THE POPULATION HAS CONVERGED TO A PRESPECIFIED SMALL PARAMETER SPACE")
  end

  # Begin evolution loops:
  lpos = 1
  criter = []
  criter_change = 1e+5

  while icall .< maxn && gnrng > x_reltol && criter_change .> f_reltol
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
            lpos = 1 + floor(npg + 0.5 - sqrt((npg + 0.5)^2 - npg * (npg + 1) * mrand()))
            idx = findall(lcs[1:k3-1] .== lpos)
            isempty(idx) && break
          end
          lcs[k3] = lpos
        end
        lcs = sort(lcs)
        # Construct the simplex:
        s = zeros(FT, nps, nopt)
        s = cx[lcs, :]
        sf = cf[lcs]

        snew, fnew, icall = cceua(fn, s, sf, bl, bu, icall)

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

    @printf("Iteration = %3d, nEvals = %3d, Best Cost = %.5f\n", nloop, icall, bestf)

    # Check for convergency
    if icall >= maxn
      exitflag = 0
      if verbose
        disp("\n*** Optimization search terminated because the limit ***")
        println("On the maximum number of trials $(maxn) has been exceeded!")
      end
    end
    if gnrng .< x_reltol
      exitflag = 1
      verbose && disp("The population has converged to a prespecified small parameter space")
    end
    push!(criter, bestf)
    # criter = [criter bestf]'
    if (nloop >= kstop)
      criter_change = abs(criter[nloop] - criter[nloop-kstop+1])
      criter_change = criter_change / mean(abs.(criter[nloop-kstop+1:nloop]))

      if criter_change .< f_reltol
        exitflag = 1
        if verbose
          println("The best point has improved in last $(num2str(kstop)) loops by less than the threshold $(num2str(f_reltol))")
          println("Convergency has achieved based on objective function criteria!!!")
        end
      end
    end
    # End of the Outer Loops
  end

  if verbose
    @printf("Search was stopped at trial number: %d \n", icall)
    println("Normalized geometric range = $(num2str(gnrng))")
    println("The best point has improved in last $(num2str(kstop)) LOOPS BY $(num2str(criter_change))")
  end
  bestx, bestf, exitflag
end
