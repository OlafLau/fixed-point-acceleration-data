using Pkg
Pkg.activate("/home/ljy/Research/Env/DW_ACSCFT")

using Polymer
using Scattering
using Polyorder
using Polyorder: SD, Anderson, SIS, ETD, NGMRES, OACCEL
import Polyorder as JP
using Makie
using CairoMakie
using DataFrames
using Random: seed!, Xoshiro, default_rng
using Statistics
using JLD2

cylinder() = AB_system(χN=18.0, fA=0.30)

function build_scft(updater)
    """
    test passed
    """
    scft = let
        ds = 0.01
        uc = UnitCell(Square(), 3.6434)
        lat = BravaisLattice(uc)
        system = cylinder()
        NoncyclicChainSCFT(system, lat, ds, updater=updater)
    end
    return scft
end

function run_single_scft!(stat_result::AbstractDict, scft::AbstractSCFT; rng=default_rng(),
                        max_iter::Int=1000, precompile::Bool=false)
    """
    test passed
    """
    if precompile
        scft_pc = deepcopy(scft)
        scft_pc.updater.warmup = 1
        scftconfig = SCFTConfig(max_iter=2)
        ioconfig = IOConfig(
            verbosity=0,
            save_config=false,
            save_w=false,
            save_ϕ=false,
            save_summary=false,
            save_trace=false
            )
        config = JP.Config(scft=scftconfig, io=ioconfig)
        JP.solve!(scft_pc, config)

        return nothing
    else
        JP.initialize!(scft, :randn, rng=rng)
        JP.reset!(scft.updater)       
        convergence = let
            scftconfig = SCFTConfig(max_iter=max_iter, tol=2*1e-5)
            ioconfig = IOConfig(
                verbosity=0,
                save_config=false,
                save_w=false,
                save_ϕ=false,
                save_summary=false,
                save_trace=false
            )
            config = JP.Config(scft=scftconfig, io=ioconfig)
            JP.solve!(scft, config)
        end
    
        # if JP.good(convergence)
        if convergence == JP.Successful()
            if isapprox(JP.F(scft), 2.587; atol=0.02)
                push!(stat_result["converge"], scft.updater.evals[end])
            else
                stat_result["converge_wrong_target"] += 1
            end
        else
            stat_result["diverge"] += 1
        end
        
        return nothing
    end
end

function create_stat_result()
    stat_result = Dict(
    "converge" => Int64[],
    "diverge" => 0,
    "converge_wrong_target" => 0,
    "summary" => Dict(
        "success_rate " => 0.0,
        "fastest" => 0,
        "slowest" => 0,
        "average" => 0.0
        )
    )

    return stat_result
end

function run_batch_stat(updater; times::Int=100)
    stat_result = create_stat_result()
    scft = build_scft(updater)
    # precompile
    run_single_scft!(stat_result, scft; precompile=true)
    # stat part : 100 times
    seed = 3083
    rng = Xoshiro(seed)
    for i in 1:times
        scft = build_scft(updater)
        run_single_scft!(stat_result, scft; rng=rng)
    end
    # summarize
    summary = stat_result["summary"]
    summary["success_rate"] = length(stat_result["converge"]) / times
    if !isempty(stat_result["converge"])
        summary["fastest"] = minimum(stat_result["converge"])
        summary["slowest"] = maximum(stat_result["converge"])
        summary["average"] = mean(stat_result["converge"])
    end

    return stat_result
end