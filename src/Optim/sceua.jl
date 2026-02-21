function sceua(fn::Function, x0::Vector{FT}, bl::Vector{FT}, bu::Vector{FT}, args...;
    verbose=true, parallel=true,
    maxn=1000, kstop=5, f_reltol=0.0001, x_reltol=0.0001,
    n_param=length(x0), # ignore this
    # 核心参数
    n_complex=5,
    size_complex=2 * n_param + 1,
    size_simplex=n_param + 1,
    n_evolu=size_complex,
    n_pop=size_complex * n_complex,
    # 其他参数
    seed=1, include_initial=1, kw...) where {FT<:Real}

    exitflag = ReturnCode.Default
    set_seed(seed)

    # Create an initial population to fill array x[n_pop, n_param]:
    bound = bu - bl

    main_rng = _rng_from_seed(_rng_seed(seed, 0, 0))
    x = zeros(FT, n_pop, n_param)
    for i = 1:n_pop
        x[i, :] = bl + rand(main_rng, FT, n_param) .* bound
    end
    include_initial == 1 && (x[1, :] = x0)

    _fn(i, _args...; _kw...) = fn(x[i, :], _args...; _kw...) |> sanitize
    xf = FT.(par_map(_fn, 1:n_pop, args...; kw..., parallel, use_deepcopy=true))

    # 判定标准
    num_evals = n_pop
    gnrng = geometric_range(x, bound)
    criter_change = 1e+5

    xf, idx = SORT(xf)
    x = x[idx, :]

    nloop = 0
    bestx = x[1, :]
    bestf = xf[1]
    @printf("Iteration = %3d, nEvals = %3d, Best Cost = %.5f\n", nloop, num_evals, bestf)

    # Begin evolution loops:
    criter = []

    eval_counts = zeros(Int, n_complex)
    nslots = parallel ? Threads.maxthreadid() : 1
    local_args = [deepcopy(args) for _ in 1:nslots]
    local_kw = [deepcopy(kw) for _ in 1:nslots]

    while num_evals .< maxn && gnrng > x_reltol && criter_change .> f_reltol
        nloop = nloop + 1

        # Loop on complexes [sub-populations]
        @par parallel for igs = 1:n_complex
            # Partition the population into complexes [sub-populations]
            k1 = 1:size_complex
            k2 = (k1 .- 1) .* n_complex .+ igs
            cx = x[k2, :]
            cf = xf[k2]

            tid = Threads.threadid()
            eval_args = local_args[tid]
            eval_kw = local_kw[tid]

            # Evolve sub-population igs for nspl steps:
            lcs = zeros(Int, size_simplex)
            local_num_evals = 0

            for loop = 1:n_evolu
                rng = _rng_from_seed(_rng_seed(seed, nloop, igs, loop))

                # Select simplex by sampling the complex according to a linear
                # probability distribution
                lcs[1] = 1
                for k3 = 2:size_simplex
                    lpos = 1
                    for iter = 1:1000
                        lpos = 1 + floor(size_complex + 0.5 - sqrt((size_complex + 0.5)^2 - size_complex * (size_complex + 1) * rand(rng)))
                        idx = findall(lcs[1:k3-1] .== lpos)
                        isempty(idx) && break
                    end
                    lcs[k3] = lpos
                end
                lcs = sort(lcs)
                # Construct the simplex:
                s = zeros(FT, size_simplex, n_param)
                s = cx[lcs, :]
                sf = cf[lcs]

                snew, fnew, local_num_evals = cceua(fn, s, sf, bl, bu, local_num_evals, eval_args...; eval_kw..., rng)

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
            eval_counts[igs] = local_num_evals
            # End of Loop on Complex Evolution
        end

        num_evals += sum(eval_counts)
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
