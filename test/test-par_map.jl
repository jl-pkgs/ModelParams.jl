@testset "par_map" begin
    A = collect(1:20)

    r1 = par_map(x -> x^2, A; parallel=false, use_deepcopy=false)
    @test r1 == map(x -> x^2, A)

    r2 = par_map(x -> x + 1, A; parallel=true, use_deepcopy=false)
    @test r2 == map(x -> x + 1, A)

    mutable = [0]
    r3 = par_map(A, mutable; parallel=false, use_deepcopy=true) do x, m
        m[1] += 1
        m[1]
    end
    @test mutable[1] == 0
    @test r3 == collect(1:length(A))

    mutable2 = [0]
    r4 = par_map(A, mutable2; parallel=true, use_deepcopy=true) do x, m
        m[1] += 1
        x
    end
    @test mutable2[1] == 0
    @test r4 == A
end



mutable struct SolverState
    counter::Int
    scratch::Vector{Float64}
end

mutable struct SolverConfig
    state::SolverState
end

function simulate_step(x, cfg::SolverConfig; delay=0.0)
    cfg.state.counter += 1
    cfg.state.scratch[1] += x
    delay > 0 && sleep(delay)
    return x * x
end



@testset "par_map thread safety and timing" begin
    A = collect(1:48)
    delay = 0.05

    cfg_safe = SolverConfig(SolverState(0, [0.0]))
    t_serial = @elapsed r_serial = par_map(simulate_step, A, cfg_safe; parallel=false, use_deepcopy=true, delay)
    @test r_serial == map(x -> x * x, A)
    @test cfg_safe.state.counter == 0
    @test cfg_safe.state.scratch[1] == 0.0

    cfg_parallel = SolverConfig(SolverState(0, [0.0]))
    t_parallel = @elapsed r_parallel = par_map(simulate_step, A, cfg_parallel; parallel=true, use_deepcopy=true)
    @test r_parallel == map(x -> x * x, A)
    @test cfg_parallel.state.counter == 0
    @test cfg_parallel.state.scratch[1] == 0.0

    speedup = t_serial / t_parallel
    println("par_map benchmark: threads=$(Threads.nthreads()), serial=$(round(t_serial, digits=4))s, parallel=$(round(t_parallel, digits=4))s, speedup=$(round(speedup, digits=2))x")

    if Threads.nthreads() > 1
        @test speedup > 1
    end
end
