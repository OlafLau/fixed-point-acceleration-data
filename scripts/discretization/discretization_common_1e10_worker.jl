include(joinpath(@__DIR__, "si_benchmark_worker.jl"))

function common_discretization_main()
    jobs = JSON.parsefile(joinpath(RUNROOT, WORKER * "_plan.json"))
    baselines = JSON.parsefile(joinpath(RUNROOT, "baseline_plan.json"))
    phases = unique(job["phase"] for job in jobs)
    algorithms = unique(job["algorithm"] for job in jobs)
    @assert length(phases) == 1
    phase = only(phases)
    warmups = Any[]

    # Exercise each updater on the inexpensive BCC system before full warm-up.
    for algorithm in algorithms
        job = copy(first(j for j in baselines if j["phase"] == "bcc" && j["algorithm"] == algorithm))
        job["id"] = "compile_" * algorithm
        job["budget"] = 6
        push!(warmups, execute(job; save=false))
    end

    # One complete excluded warm-up for every measured algorithm on this phase.
    for source in baselines
        source["phase"] == phase || continue
        job = merge(copy(source), Dict("id" => "warmup_" * source["id"]))
        atomicjson(joinpath(RUNROOT, "workers", WORKER * ".json"),
            Dict("state" => "warming", "id" => job["id"], "pid" => getpid()))
        row = execute(job; save=false)
        push!(warmups, row)
        atomicjson(joinpath(RUNROOT, "warmups", WORKER * ".json"), warmups)
        @assert row["target_validated"] "Baseline candidate failed target validation"
        println("FULL_WARMUP ", job["id"], " N=", row["actual_neval"])
        flush(stdout)
    end

    for job in jobs
        atomicjson(joinpath(RUNROOT, "workers", WORKER * ".json"),
            Dict("state" => "measuring", "id" => job["id"], "pid" => getpid()))
        row = execute(job)
        println("RESULT ", job["id"], " ", row["status"], " N=", row["actual_neval"],
            " T=", row["solve_seconds"], " eligible=", row["timing_eligible"])
        flush(stdout)
    end
    atomicjson(joinpath(RUNROOT, "workers", WORKER * ".json"),
        Dict("state" => "complete", "pid" => getpid()))
end

common_discretization_main()
