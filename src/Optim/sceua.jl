"""
    sceua(fn::Function, x0::Vector, bl::Vector, bu::Vector;
        maxn=500, kstop=5, pcento=0.01, peps=0.001, ngs=5, iseed=-1, iniflg=1)

# Parameters:
- `x0`       : the initial parameter array at the start(), the optimized parameter array at the end
- `fn`       : the objective function value corresponding to the optimized parameters
- `bl`       : the lower bound of the parameters
- `bu`       : the upper bound of the parameters
- `iseed`    : the random seed number [for repetetive testing purpose]
- `iniflg`   : flag for initial parameter array (default 1, included it in initial population; otherwise, not included)
- `ngs`      : number of complexes [sub-populations]
- `npg`      : number of members in a complex()
- `nps`      : number of members in a simplex
- `nspl`     : number of evolution steps for each complex before shuffling
- `mings`    : minimum number of complexes required during the optimization process
- `maxn`     : maximum number of function evaluations allowed during optimization
- `kstop`    : maximum number of evolution loops before convergency
- `percento` : the percentage change allowed in kstop loops before convergency
  
# Return:
`bestx, bestf, exitflag, output`

# Examples:
```julia
x, feval, exitflag, output = sceua(fn, x0, bl, bu)
```
"""
function sceua(fn::Function, x0::Vector{FT}, bl::Vector{FT}, bu::Vector{FT};
  verbose=false,
  maxn=1000, kstop=5, pcento=0.01, peps=0.0001, ngs=5, iseed=1, iniflg=1) where {FT<:Real}

  exitflag = -1
  set_seed(iseed)

  # Initialize SCE parameters:
  nopt = length(x0)
  npg = 2 * nopt + 1
  nps = nopt + 1
  nspl = npg
  npt = npg * ngs

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

  # worstx = x[npt, :]
  # worstf = xf[npt]
  # BESTF = bestf
  # BESTX = bestx
  # ICALL = icall

  # Compute the standard deviation for each parameter
  xnstd = std(x)
  # Computes the normalized geometric range of the parameters
  ## TODO: BUG HERE
  gnrng = exp(mean(log.((colMax(x) - colMin(x)) ./ bound)))

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

  if gnrng .< peps
    exitflag = 1
    disp("THE POPULATION HAS CONVERGED TO A PRESPECIFIED SMALL PARAMETER SPACE")
  end

  # Begin evolution loops:
  lpos = 1
  criter = []
  criter_change = 1e+5

  while icall .< maxn && gnrng > peps && criter_change .> pcento
    nloop = nloop + 1
    # Loop on complexes [sub-populations]
    for igs = 1:ngs
      # Partition the population into complexes [sub-populations]
      k1 = 1:npg
      k2 = (k1 .- 1) .* ngs .+ igs
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
            if isempty(idx)
              break
            end
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

    # Compute the standard deviation for each parameter
    xnstd = std(x)
    # Computes the normalized geometric range of the parameters
    gnrng = exp(mean(log.((colMax(x) - colMin(x)) ./ bound)))

    # disp(["Evolution Loop: ' num2str(nloop) ', Trial: " num2str(icall)])
    # disp(["BESTF  : ' num2str(bestf), ' ' 'BESTX  : [' num2str(bestx), ']"])
    # disp(["WORSTF : ' num2str(worstf), ' ' 'WORSTX : [' num2str(worstx), ']"])
    # disp(' ')
    @printf("Iteration = %3d, nEvals = %3d, Best Cost = %.5f\n", nloop, icall, bestf)

    # Check for convergency
    if icall >= maxn
      exitflag = 0
      if verbose
        disp("\n*** Optimization search terminated because the limit ***")
        println("On the maximum number of trials $(maxn) has been exceeded!")
      end
    end
    if gnrng .< peps
      exitflag = 1
      verbose && disp("The population has converged to a prespecified small parameter space")
    end
    push!(criter, bestf)
    # criter = [criter bestf]'
    if (nloop >= kstop)
      criter_change = abs(criter[nloop] - criter[nloop-kstop+1]) * 100
      criter_change = criter_change / mean(abs.(criter[nloop-kstop+1:nloop]))

      if criter_change .< pcento
        exitflag = 1
        if verbose
          println("The best point has improved in last $(num2str(kstop)) loops by less than the threshold $(num2str(pcento))")
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
