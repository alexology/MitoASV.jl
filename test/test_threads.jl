@testset "Threaded loops wait for every task" begin
    if Threads.nthreads() == 1
        @test_skip "needs more than one thread"
    else
        # Index 1 fails at once, index 2 blocks until released: the caller
        # must not return (or throw) while the blocked task is still running.
        function run_blocked(call)
            release = Channel{Nothing}(1)
            finished = Threads.Atomic{Bool}(false)
            block() = (take!(release); finished[] = true)
            caller = Threads.@spawn try
                call(block)
                nothing
            catch e
                e
            end
            sleep(0.5)
            early = istaskdone(caller)
            put!(release, nothing)
            err = fetch(caller)
            (; early, finished = finished[], err)
        end
        failing(i) = error("failure at $i")

        r = run_blocked(block -> MitoASV._tforeach(3) do i
            i == 1 && failing(i)
            i == 2 && block()
        end)
        @test !r.early
        @test r.finished
        @test r.err isa TaskFailedException
        @test occursin("failure at 1", sprint(showerror, r.err))

        r = run_blocked(block -> MitoASV.threshold_edges(3, 1.0) do i, j, t
            i == 1 && failing(i)
            i == 2 && block()
            0.0
        end)
        @test !r.early
        @test r.finished
        @test r.err isa TaskFailedException
        @test occursin("failure at 1", sprint(showerror, r.err))

        # Without failures the results are unchanged.
        out = zeros(Int, 50)
        MitoASV._tforeach(i -> (out[i] = i^2), 50)
        @test out == (1:50) .^ 2
        @test MitoASV.threshold_edges((i, j, t) -> abs(i - j) / 10, 6, 0.2) ==
              [(i, j, (j - i) / 10) for i in 1:6 for j in (i + 1):min(6, i + 2)]
    end
end
