# BCC time-per-evaluation benchmark vs history length m, for the nine accelerators.
# Reuses the principal-campaign BCC setup and instrumentation (process-local), scanning m
# at each method's BCC-optimal step size. Usage: julia timing_m_bcc.jl <shard> <nshards> [reps]
using Polyorder, Polymer, FFTW, LinearAlgebra, Statistics, CSV, DataFrames, HDF5
import Polyorder as JP
@assert pkgversion(JP)==v"0.26.7"

const OUT = get(ENV, "SCFT_OUTPUT_DIR", joinpath(@__DIR__, "results")); mkpath(OUT)
const INPUT = get(ENV, "SCFT_INPUT_DIR", joinpath(@__DIR__, "inputs"))
FFTW.set_num_threads(1); BLAS.set_num_threads(1)

const NQ = Ref(0); const QTIME = Ref(0.0); const LIMIT = Ref(6000)
struct BudgetExhausted <: Exception end
@eval JP function q!(scft::NoncyclicChainSCFT)
    Main.NQ[] >= Main.LIMIT[] && throw(Main.BudgetExhausted())
    t = time_ns()
    update_propagator!(scft); update_density!(scft)
    q!(scft.forces, scft.wfields, scft.ϕfields, scft.system.χNmatrix)
    Main.NQ[] += 1; Main.QTIME[] += (time_ns()-t)/1e9
    return nothing
end

function ab3(fa, chi)
    a = KuhnSegment(:A); b = KuhnSegment(:B); junction = BranchPoint(:EB)
    blocks = [PolymerBlock(:A,a,fa,FreeEnd(:A),junction),
              PolymerBlock(:B1,b,(1-fa)/3,junction,FreeEnd(:B1)),
              PolymerBlock(:B2,b,(1-fa)/3,junction,FreeEnd(:B2)),
              PolymerBlock(:B3,b,(1-fa)/3,junction,FreeEnd(:B2))]
    return PolymerSystem([Component(BlockCopolymer(:AB3,blocks))], Dict([:A,:B]=>chi))
end

function updater(name, alpha, m)
    base = endswith(name,"etd") ? JP.ETD(alpha) : endswith(name,"sis") ? JP.SIS(alpha) : JP.SD(alpha)
    startswith(name,"aa_") && return JP.Anderson(base; m, warmup=0, αw=.2)
    startswith(name,"na_") && return JP.NGMRES(base;  m, warmup=0, αw=.2)
    startswith(name,"oa_") && return JP.OACCEL(base;  m, warmup=0, αw=.2)
    return base
end

function build(alg, alpha, m)
    lat = BravaisLattice(UnitCell(Cubic(), 3.900))
    w = JP.AuxiliaryField(zeros(27,27,27), lat)
    sc = NoncyclicChainSCFT(ab3(0.4, 16.0), w, 0.01; mde=OSF,
                            updater=updater(alg,alpha,m), init=:zeros)
    JP.initialize!(sc, joinpath(INPUT, "bcc_seed.h5"))
    return sc
end

# BCC-optimal step size (from the principal campaign, neval objective)
const SPEC = Dict("aa_etd"=>1.0,"aa_sd"=>0.8,"aa_sis"=>1.0,
                  "na_etd"=>1.0,"na_sd"=>0.2,"na_sis"=>1.0,
                  "oa_etd"=>1.0,"oa_sd"=>0.2,"oa_sis"=>0.6)
const MS = (10,20,30,50,80,100)

function runone(alg, m, rep; budget=6000)
    sc = build(alg, SPEC[alg], m)
    cfg = JP.Config(scft=SCFTConfig(max_iter=budget, tol=1e-10, norm=:norm2, norm2=:mean,
                                    relative=true, maxtime=1e8),
                    io=IOConfig(verbosity=0, display_interval=1000000, save_config=false,
                                save_w=false, save_ϕ=false, save_summary=false, save_trace=false))
    GC.gc(); NQ[]=0; QTIME[]=0.0; LIMIT[]=budget
    t0 = time_ns(); status = "error"
    try
        result = JP.solve!(sc, cfg)
        status = result == JP.Successful() ? "converged" : "solver_stopped"
    catch e
        status = e isa BudgetExhausted ? "budget_exhausted" : "numerical_breakdown"
    end
    sec = (time_ns()-t0)/1e9
    r = isempty(sc.updater.rs) ? Inf : sc.updater.rs[end]
    return (algorithm=alg, alpha=SPEC[alg], m=m, rep=rep, status=status,
            neval=NQ[], seconds=sec, q_seconds=QTIME[], residual=r)
end

function main()
    shard = parse(Int, ARGS[1]); nsh = parse(Int, ARGS[2])
    reps = length(ARGS)>2 ? parse(Int,ARGS[3]) : 5
    0 <= shard < nsh || error("shard must satisfy 0 <= shard < nsh")
    # Interleave complete parameter sweeps by repetition.  Assign repetitions,
    # rather than parameter combinations, to shards so every (algorithm, m)
    # pair is sampled equally on every NUMA node.
    jobs = [(a,m,r) for r in 1:reps for a in sort(collect(keys(SPEC))) for m in MS]
    mine = [j for j in jobs if (j[3]-1) % nsh == shard]
    out = joinpath(OUT, "bcc_timing_results$(nsh>1 ? "_shard"*string(shard) : "").csv")
    # A complete, unrecorded solve for each algorithm compiles the solver,
    # accelerator, preconditioner, and history-management paths before timing.
    # In particular, this avoids counting Julia's first-call JIT cost in rep=1.
    for a in unique(first.(mine))
        warm = runone(a, maximum(MS), 0)
        warm.status == "converged" || error("JIT warm-up failed for $a: $(warm.status)")
    end
    rows = NamedTuple[]
    for (a,m,r) in mine
        push!(rows, runone(a,m,r))
        length(rows) % 10 == 0 && (CSV.write(out, DataFrame(rows)); println("shard",shard," ",length(rows),"/",length(mine)); flush(stdout))
    end
    CSV.write(out, DataFrame(rows))
    println("saved ",length(rows)," to ",out)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
