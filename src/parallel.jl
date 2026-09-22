export @par, par_map, get_clusters, get_nthreads

using Base.Threads
using ProgressMeter

macro par(parallel, ex)
    ex isa Expr && ex.head == :for ||
        throw(ArgumentError("@par requires a for-loop expression"))

    esc(quote
        if $parallel
            @threads $ex
        else
            $ex
        end
    end)
end

get_clusters() = Threads.nthreads()
get_nthreads() = Threads.nthreads(:interactive) + Threads.nthreads(:default)

function par_map(fn, A, args...;
    parallel::Bool=true, progress::Bool=false, use_deepcopy::Bool=true, kw...)

    res = Array{Any}(undef, size(A))

    total = length(A)
    nworkers = parallel ? get_clusters() : 1
    chunk_size = max(1, cld(total, nworkers))
    chunks = collect(Iterators.partition(eachindex(A), chunk_size))

    p = progress ? Progress(total) : nothing
    @par parallel for chunk in chunks
        eval_args = use_deepcopy ? deepcopy(args) : args
        eval_kw = use_deepcopy ? deepcopy(kw) : kw

        for i in chunk
            res[i] = fn(A[i], eval_args...; eval_kw...)
            progress && next!(p)
        end
    end
    map(x -> x, res) # unpack, which necessary
end
