using Polyorder, Polymer, Scattering, FFTW, LinearAlgebra, Statistics, CSV, DataFrames, JLD2, JSON, SHA, HDF5, Logging
import Polyorder as JP
@assert pkgversion(JP)==v"0.26.7"
const RUNROOT=abspath(ARGS[1])
const WORKER=length(ARGS)>1 ? ARGS[2] : "serial"
FFTW.set_num_threads(1); BLAS.set_num_threads(1)
global_logger(ConsoleLogger(stderr,Logging.Warn))
const INPUT=joinpath(RUNROOT,"inputs")
const NQ=Ref(0); const QTIME=Ref(0.0); const LIMIT=Ref(6000)
struct BudgetExhausted <: Exception end
@eval JP function q!(scft::NoncyclicChainSCFT)
    Main.NQ[] >= Main.LIMIT[] && throw(Main.BudgetExhausted())
    t=time_ns()
    update_propagator!(scft)
    update_density!(scft)
    q!(scft.forces,scft.wfields,scft.ϕfields,scft.system.χNmatrix)
    Main.NQ[]+=1;Main.QTIME[]+=(time_ns()-t)/1e9
    return nothing
end

function ab3(fa,chi)
    a=KuhnSegment(:A);b=KuhnSegment(:B);junction=BranchPoint(:EB)
    # Preserve the recovered factory exactly, including its terminal labels.
    blocks=[PolymerBlock(:A,a,fa,FreeEnd(:A),junction),
        PolymerBlock(:B1,b,(1-fa)/3,junction,FreeEnd(:B1)),
        PolymerBlock(:B2,b,(1-fa)/3,junction,FreeEnd(:B2)),
        PolymerBlock(:B3,b,(1-fa)/3,junction,FreeEnd(:B2))]
    return PolymerSystem([Component(BlockCopolymer(:AB3,blocks))],Dict([:A,:B]=>chi))
end

function updater(name,alpha,m)
    base=endswith(name,"etd") ? JP.ETD(alpha) : endswith(name,"sis") ? JP.SIS(alpha) : JP.SD(alpha)
    startswith(name,"aa_") && return JP.Anderson(base;m,warmup=0,αw=.2)
    startswith(name,"na_") && return JP.NGMRES(base;m,warmup=0,αw=.2)
    startswith(name,"oa_") && return JP.OACCEL(base;m,warmup=0,αw=.2)
    return base
end

function build(job)
    phase=job["phase"]; fa=phase=="bcc" ? .4 : .34;chi=phase=="bcc" ? 16. : 20.
    cell=phase=="bcc" ? UnitCell(Cubic(),3.900) : UnitCell(Tetragonal(),12.3042,6.501)
    lat=BravaisLattice(cell)
    grid=Int.(get(job,"grid",phase=="bcc" ? [27,27,27] : [84,84,45]))
    w=JP.AuxiliaryField(zeros(grid...),lat)
    sc=NoncyclicChainSCFT(ab3(fa,chi),w,get(job,"ds",.01);mde=OSF,
        updater=updater(job["algorithm"],job["alpha"],job["m"]),init=:zeros)
    JP.initialize!(sc,joinpath(INPUT,phase*"_seed.h5"))
    return sc
end

const REFS=Dict{String,Array{Float64,3}}()
function reference(phase,shape)
    key=phase*join(shape,"x")
    get!(REFS,key) do
        rho=h5read(joinpath(INPUT,phase*"_reference_density.h5"),"phiA")
        size(rho)==shape ? rho : JP.resample(rho,shape)
    end
end

function descriptors(sc,phase)
    rho=sc.ϕfields[1].data; ref=reference(phase,size(rho))
    a=rho .- mean(rho); b=ref .- mean(ref)
    fa=fft(a);fb=fft(b); pa=abs.(fa);pb=abs.(fb)
    spectral=dot(vec(pa),vec(pb))/max(norm(pa)*norm(pb),eps())
    # Translation-registered correlation and normalized RMS error; no rotation fit.
    cross=maximum(real.(ifft(fa.*conj.(fb))))
    nrmse=sqrt(max(sum(abs2,a)+sum(abs2,b)-2cross,0.))/max(norm(b),eps())
    return std(rho),mean(rho),spectral,nrmse
end

function atomicjson(path,value)
    temp=path*".tmp."*WORKER
    open(temp,"w") do io;JSON.print(io,value,2);end
    mv(temp,path;force=true)
end

function runjob(job;save=true)
    id=job["id"]; phase=job["phase"]
    started=time_ns(); sc=build(job);setup=(time_ns()-started)/1e9
    cfg=JP.Config(scft=SCFTConfig(max_iter=job["budget"],tol=get(job,"tol",1e-10),norm=:norm2,
        norm2=:mean,relative=true,maxtime=1e8),
        io=IOConfig(verbosity=0,display_interval=1000000,save_config=false,save_w=false,
            save_ϕ=false,save_summary=false,save_trace=false))
    GC.gc();NQ[]=0;QTIME[]=0;LIMIT[]=job["budget"]
    solve_start=time_ns();status="error";msg="";result=nothing
    try
        result=JP.solve!(sc,cfg)
        status=result==JP.Successful() ? "converged" : "solver_stopped"
    catch e
        status=e isa BudgetExhausted ? "budget_exhausted" : "numerical_breakdown"
        msg=sprint(showerror,e)
    end
    seconds=(time_ns()-solve_start)/1e9
    r=isempty(sc.updater.rs) ? Inf : sc.updater.rs[end]
    energy=isempty(sc.updater.Fs) ? NaN : sc.updater.Fs[end]
    rms=NaN;inc=NaN;amp=NaN;meanrho=NaN;similarity=NaN;nrmse=NaN
    accepted=status=="converged" && isfinite(r) && r<get(job,"tol",1e-10)
    if accepted
        ratios=[norm(sc.forces[i])/max(norm(sc.wfields[i]),eps()) for i in 1:2]
        rms=sqrt(mean(abs2,ratios));inc=maximum(abs,sc.ϕfields[1].data+sc.ϕfields[2].data .-1)
        amp,meanrho,similarity,nrmse=descriptors(sc,phase)
    end
    refF=phase=="bcc" ? 2.8338924026372423 : 3.397773069273983
    target=accepted && abs(energy-refF)<2e-6 && similarity>.99 && nrmse<.05
    row=merge(copy(job),Dict("status"=>status,"solver_status"=>string(result),"message"=>msg,
        "worker"=>WORKER,"polyorder"=>string(pkgversion(JP)),"julia"=>string(VERSION),
        "julia_threads"=>Threads.nthreads(),"fft_threads"=>FFTW.get_num_threads(),"blas_threads"=>BLAS.get_num_threads(),
        "actual_neval"=>NQ[],"stored_neval"=>isempty(sc.updater.evals) ? 0 : sc.updater.evals[end],
        "iterations"=>length(sc.updater.rs),"free_energy"=>energy,"residual_mean"=>r,"residual_rms"=>rms,
        "incompressibility"=>inc,"density_sd"=>amp,"density_mean"=>meanrho,"spectral_cosine"=>similarity,
        "registered_nrmse"=>nrmse,"target_validated"=>target,"solve_seconds"=>seconds,
        "q_seconds"=>QTIME[],"setup_seconds"=>setup,"grid"=>collect(size(sc.wfields[1])),"ds"=>get(job,"ds",.01)))
    if save
        CSV.write(joinpath(RUNROOT,"traces",id*".csv"),DataFrame(update=1:length(sc.updater.rs),
            stored_neval=sc.updater.evals,residual=sc.updater.rs,free_energy=sc.updater.Fs))
        if accepted && get(job,"savefields",false)
            jldsave(joinpath(RUNROOT,"fields",id*".jld2");w=[copy(x.data) for x in sc.wfields],
                density=[copy(x.data) for x in sc.ϕfields],metadata=row)
        end
        atomicjson(joinpath(RUNROOT,"results",id*".json"),row)
    end
    return row
end

function main()
    jobs=JSON.parsefile(joinpath(RUNROOT,"plan.json"))
    names=unique(j["algorithm"] for j in jobs)
    atomicjson(joinpath(RUNROOT,"workers",WORKER*".json"),Dict("state"=>"compiling","pid"=>getpid()))
    for name in names
        runjob(Dict("id"=>"compile","phase"=>"bcc","algorithm"=>name,"alpha"=>.4,"m"=>10,"budget"=>6);save=false)
    end
    for job in jobs
        id=job["id"]
        isfile(joinpath(RUNROOT,"results",id*".json")) && continue
        claim=joinpath(RUNROOT,"claims",id)
        try;mkdir(claim);catch e;isdir(claim) ? continue : rethrow();end
        atomicjson(joinpath(RUNROOT,"workers",WORKER*".json"),Dict("state"=>"running","id"=>id,"pid"=>getpid()))
        row=runjob(job)
        println(id," ",row["status"]," N=",row["actual_neval"]," target=",row["target_validated"]);flush(stdout)
    end
    atomicjson(joinpath(RUNROOT,"workers",WORKER*".json"),Dict("state"=>"finished","pid"=>getpid()))
end
if abspath(PROGRAM_FILE)==@__FILE__;main();end
