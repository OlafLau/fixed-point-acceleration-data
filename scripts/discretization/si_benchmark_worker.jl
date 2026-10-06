include(joinpath(@__DIR__, "principal_worker.jl"))
include(joinpath(@__DIR__, "discretization_registration.jl"))
using Dates
@assert VERSION==v"1.11.5"
const CPU=parse(Int,ARGS[3]);const NODE=parse(Int,ARGS[4])
const DIRECT=Int[]
const REFINFO=JSON.parsefile(joinpath(RUNROOT,"references.json"))
utcstamp()=Dates.format(now(UTC),dateformat"yyyy-mm-ddTHH:MM:SS.sss")*"Z"
@eval JP function train!(model::SCFTModel,n)
    foreach(1:n) do _
        ProgressMeter.next!(model.prog)
        model.config.scft.symmetrize && symmetrize!(model.scft)
        update!(model.scft,model.updater,model.config)
        push!(Main.DIRECT,Main.NQ[])
    end
    return nothing
end

function execute(job;save=true)
    started=utcstamp();t=time_ns();sc=build(job);setup=(time_ns()-t)/1e9
    cfg=JP.Config(scft=SCFTConfig(max_iter=job["budget"],tol=job["tol"],norm=:norm2,norm2=:mean,
        relative=true,maxtime=1e8),io=IOConfig(verbosity=0,display_interval=1000000,
        save_config=false,save_w=false,save_ϕ=false,save_summary=false,save_trace=false))
    contours=[]
    for (ic,qs) in enumerate(sc.propagators),(direction,q) in qs
        @assert isapprox(q.ds*(q.Ns-1),q.block.f;atol=1e-14)
        push!(contours,Dict("component"=>ic,"direction"=>string(direction),"block"=>string(Polymer.label(q.block)),
            "f"=>q.block.f,"actual_ds"=>q.ds,"nodes"=>q.Ns))
    end
    grid=collect(size(sc.wfields[1]));edges=collect(sc.wfields[1].lattice.unitcell.edges)
    @assert grid==job["grid"]
    @assert edges ≈ (job["phase"]=="bcc" ? [3.9,3.9,3.9] : [12.3042,12.3042,6.501])
    affinity=strip(match(r"(?m)^Cpus_allowed_list:\s*(.+)$",read("/proc/self/status",String)).captures[1])
    @assert affinity==string(CPU)
    GC.gc();NQ[]=0;QTIME[]=0;LIMIT[]=job["budget"];empty!(DIRECT)
    Base.cumulative_compile_timing(true);before=Base.cumulative_compile_time_ns();t=time_ns()
    result=nothing;status="error";message=""
    try
        result=JP.solve!(sc,cfg)
        status=result==JP.Successful() ? "converged" : "solver_stopped"
    catch e
        status=e isa BudgetExhausted ? "budget_exhausted" : "numerical_breakdown"
        message=sprint(showerror,e)
    end
    seconds=(time_ns()-t)/1e9;after=Base.cumulative_compile_time_ns();Base.cumulative_compile_timing(false)
    compile=(after[1]-before[1])/1e9
    expected=startswith(job["algorithm"],"aa_") ? collect(1:length(DIRECT)).+1 : 2 .* collect(1:length(DIRECT)).+2
    @assert DIRECT==expected && length(DIRECT)==length(sc.updater.rs)
    correction=startswith(job["algorithm"],"na_") ? 3 : 0
    @assert sc.updater.evals .+ correction==DIRECT
    @assert NQ[]<=job["budget"]
    r=isempty(sc.updater.rs) ? NaN : sc.updater.rs[end]
    energy=isempty(sc.updater.Fs) ? NaN : sc.updater.Fs[end]
    accepted=status=="converged" && isfinite(r) && r<job["tol"]
    nrmse=NaN;cosine=NaN;df=NaN;amp=NaN;target=false
    ref=REFINFO[job["case_id"]]
    if accepted
        @assert DIRECT[end]==NQ[]
        rho=sc.ϕfields[1].data;amp=std(rho)
        if ref["policy"]=="same_discretization"
            rid=ref["id"]
            rd=JSON.parsefile(joinpath(RUNROOT,"references",rid*".json"))
            @assert rd["grid"]==grid && rd["ds"]==job["ds"] && rd["phase"]==job["phase"]
            rr=JLD2.load(joinpath(RUNROOT,"references",rid*".jld2"),"density")[1]
            df=abs(energy-rd["free_energy"])
        else
            rr=reference(job["phase"],size(rho))
        end
        nrmse,cosine,_=integer_registered_metrics(rho,rr)
        target=amp>1e-6 && nrmse<.05 && cosine>.99 && (ref["policy"]!="same_discretization" || df<2e-6)
    end
    row=merge(copy(job),Dict("status"=>status,"solver_status"=>string(result),"message"=>message,
        "numerically_converged"=>accepted,"target_validated"=>target,"audit_passed"=>true,
        "actual_neval"=>NQ[],"iterations"=>length(DIRECT),"direct_step_trace"=>true,
        "residual_mean"=>r,"free_energy"=>energy,"registered_nrmse"=>nrmse,"spectral_cosine"=>cosine,
        "reference_energy_difference"=>df,"reference_policy"=>ref["policy"],"density_sd"=>amp,
        "grid"=>grid,"n_grid"=>prod(grid),"cell_edges"=>edges,"actual_contour"=>contours,
        "solve_seconds"=>seconds,"q_seconds"=>QTIME[],"setup_seconds"=>setup,
        "compile_seconds"=>compile,"recompile_seconds"=>(after[2]-before[2])/1e9,
        "timing_eligible"=>save && target && compile==0.0,
        "timing_role"=>save ? job["stage"] : "excluded_full_warmup",
        "julia"=>string(VERSION),"polyorder"=>string(pkgversion(JP)),
        "julia_threads"=>Threads.nthreads(),"fft_threads"=>FFTW.get_num_threads(),"blas_threads"=>BLAS.get_num_threads(),
        "cpu_affinity"=>CPU,"numa_node"=>NODE,"pid"=>getpid(),"worker"=>WORKER,
        "started_utc"=>started,"ended_utc"=>utcstamp(),"startup_convention"=>"stock_0.26.7"))
    if save
        id=job["id"]
        CSV.write(joinpath(RUNROOT,"traces",id*".csv"),DataFrame(update=1:length(DIRECT),
            actual_neval=copy(DIRECT),stored_neval=sc.updater.evals,residual=sc.updater.rs,free_energy=sc.updater.Fs))
        if accepted && get(job,"savefields",true)
            jldsave(joinpath(RUNROOT,"fields",id*".jld2");w=[copy(x.data) for x in sc.wfields],
                density=[copy(x.data) for x in sc.ϕfields],metadata=row)
        end
        atomicjson(joinpath(RUNROOT,"results",id*".json"),row)
    end
    return row
end

function main()
    jobs=JSON.parsefile(joinpath(RUNROOT,WORKER*"_plan.json"))
    baselines=JSON.parsefile(joinpath(RUNROOT,"baseline_plan.json"))
    warmups=[]
    for alg in unique(j["algorithm"] for j in jobs)
        j=copy(first(j for j in baselines if j["phase"]=="bcc" && j["algorithm"]==alg))
        j["id"]="compile_"*alg;j["budget"]=6
        push!(warmups,execute(j;save=false))
    end
    for job in baselines
        j=merge(copy(job),Dict("id"=>"warmup_"*job["id"]))
        atomicjson(joinpath(RUNROOT,"workers",WORKER*".json"),Dict("state"=>"warming","id"=>j["id"],"pid"=>getpid()))
        row=execute(j;save=false);push!(warmups,row)
        atomicjson(joinpath(RUNROOT,"warmups",WORKER*".json"),warmups)
        @assert row["target_validated"] "Baseline candidate failed target validation; stop rather than silently retune"
        println("FULL_WARMUP ",j["id"]," N=",row["actual_neval"]);flush(stdout)
    end
    for job in jobs
        atomicjson(joinpath(RUNROOT,"workers",WORKER*".json"),Dict("state"=>"measuring","id"=>job["id"],"pid"=>getpid()))
        row=execute(job)
        println("RESULT ",job["id"]," ",row["status"]," N=",row["actual_neval"]," T=",row["solve_seconds"]," eligible=",row["timing_eligible"]);flush(stdout)
    end
    atomicjson(joinpath(RUNROOT,"workers",WORKER*".json"),Dict("state"=>"complete","pid"=>getpid()))
end
if abspath(PROGRAM_FILE)==@__FILE__;main();end
