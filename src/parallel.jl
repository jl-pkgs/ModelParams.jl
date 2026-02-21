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


function par_map(fn, A, args...;
    parallel::Bool=true, progress::Bool=false, use_deepcopy::Bool=true, kw...)

    res = Array{Any}(undef, size(A))

    nslots = use_deepcopy ? Threads.maxthreadid() : 0
    local_args = use_deepcopy ? [deepcopy(args) for _ in 1:nslots] : nothing
    local_kw = use_deepcopy ? [deepcopy(kw) for _ in 1:nslots] : nothing

    total = length(A)
    p = progress ? Progress(total) : nothing
    @par parallel (
        for i in eachindex(A)
            tid = threadid()
            eval_args = use_deepcopy ? local_args[tid] : args
            eval_kw = use_deepcopy ? local_kw[tid] : kw

            res[i] = fn(A[i], eval_args...; eval_kw...)
            progress && next!(p)
        end
    )
    map(x -> x, res) # unpack, which necessary
end

get_clusters() = Threads.nthreads()


export @par, par_map, get_clusters
