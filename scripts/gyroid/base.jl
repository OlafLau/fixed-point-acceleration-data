include(joinpath(@__DIR__, "init.jl"))
SEED = get(ENV, "GYROID_REFERENCE_SEED", joinpath(@__DIR__, "inputs", "fields.h5"))

function getseed_config(xN::Float64, fA::Float64)
    ioconfig = IOConfig(
        verbosity=-1,
        save_config=false,
        save_ϕ=false,
        save_summary=false,
        save_trace=false,
        fields="seeds/xN$(xN)_fA$(fA)_fields",
        save_w=true
    )
    scftconfig = SCFTConfig(tol=1e-10)
    config = JP.Config(io=ioconfig, scft=scftconfig)
    return config
end

function fixcell_config()
    ioconfig = IOConfig(
        verbosity=-1,
        save_config=false,
        save_ϕ=false,
        save_summary=false,
        save_trace=false,
        save_w=false
    )
    scftconfig = SCFTConfig(tol=1e-8)
    config = JP.Config(io=ioconfig, scft=scftconfig)
    return config
end

function get_seed(xN::Float64, fA::Float64)
    updater = Anderson(SD(0.2); m=50)
    scft = let
        system = AB3_system(fA, xN)
        AB3_gyr(system, updater)
    end
    JP.initialize!(scft, SEED)
    config = fixcell_config()
    JP.solve!(scft, config)
    
    # cell solve
    scft_vc = JP.clone(scft);
    config = getseed_config(xN, fA)
    converged = JP.cell_solve!(scft_vc, VariableCell(BB(0.6), Anderson(SD(0.2); m=10)), config);
    accuracy = (converged[1] == JP.Successful())

    return accuracy, scft_vc.wfields
end
# for a given seedfile (.h5 format), calculate its initial residual for (fA0.5, xN30) system
function cal_init_residual(wfields)
    updater = SD(0.2)
    scft = AB3_gyr(updater)
    JP.initialize!(scft, wfields)
    ioconfig = IOConfig(verbosity=0,
        save_config=false,
        save_ϕ=false,
        save_summary=false,
        save_trace=false,
        save_w=false)
    scftconfig = SCFTConfig(max_iter=2)
    config = JP.Config(scft=scftconfig, io=ioconfig)
    JP.solve!(scft, config)

    return scft.updater.rs[1]
end
