α = [0.4, 0.6, 0.8, 1.0]
α_sd = [0.2, 0.3, 0.4, 0.5]
m = [30, 50, 80, 100]
params = [(a, m) for a in α, m in m]

function build_scft(updater::SCFTAlgorithm, seedfile::String)
    scft = AB3_gyr(updater)
    JP.initialize!(scft, seedfile)
    return scft
end

function my_ioconfig()
    ioconfig = IOConfig(
        verbosity=-1,
        save_config=false,
        save_ϕ=false,
        save_summary=false,
        save_trace=false,
        save_w=false
    )
    return ioconfig
end

function sd_stat(csv_dir::String, seedfile::String)
    sd_df = DataFrame(α=Float64[], m=Int64[], converged=Bool[], free_energy=Float64[], fevals=Int64[])
    @showprogress for a in α_sd
        updater= SD(a)
        scft = build_scft(updater, seedfile)
        ioconfig = my_ioconfig()
        scftconfig = SCFTConfig(max_iter=3000, tol=1e-8)
        config = JP.Config(scft=scftconfig, io=ioconfig)
        convergence = JP.solve!(scft, config)
        converged = (JP.Successful() == convergence)
        free_energy = JP.F(scft)
        push!(sd_df, (a, 0, converged, free_energy, scft.updater.evals[end]))
    end
    CSV.write(joinpath(csv_dir, "sd_benchmark.csv"), sd_df)
end

function anderson_sd_stat(csv_dir::String, seedfile::String)
    anderson_sd_df = DataFrame(α=Float64[], m=Int64[], converged=Bool[], free_energy=Float64[], fevals=Int64[])
    @showprogress for (a, m) in params
        updater = Anderson(SD(a); m=m, warmup=0)
        scft = build_scft(updater, seedfile)
        ioconfig = my_ioconfig()
        scftconfig = SCFTConfig(max_iter=1000, tol=1e-8)
        config = JP.Config(scft=scftconfig, io=ioconfig)
        convergence = JP.solve!(scft, config)
        converged = (JP.Successful() == convergence)
        free_energy = JP.F(scft) 
        push!(anderson_sd_df, (a, m, converged, free_energy, scft.updater.evals[end]))    
    end
    CSV.write(joinpath(csv_dir, "anderson_sd_benchmark.csv"), anderson_sd_df)
end

function oaccel_sd_stat(csv_dir::String, seedfile::String)
    oaccel_sd_df = DataFrame(α=Float64[], m=Int64[], converged=Bool[], free_energy=Float64[], fevals=Int64[])
    @showprogress for (a, m) in params
        updater = OACCEL(SD(a); m=m, warmup=0)
        scft = build_scft(updater, seedfile)
        ioconfig = my_ioconfig()
        scftconfig = SCFTConfig(max_iter=1000, tol=1e-8)
        config = JP.Config(scft=scftconfig, io=ioconfig)
        convergence = JP.solve!(scft, config)
        converged = (JP.Successful() == convergence)
        free_energy = JP.F(scft) 
        push!(oaccel_sd_df, (a, m, converged, free_energy, scft.updater.evals[end]))    
    end
    CSV.write(joinpath(csv_dir, "oaccel_sd_benchmark.csv"), oaccel_sd_df)
end

function ngmres_sd_stat(csv_dir::String, seedfile::String)
    ngmres_sd_df = DataFrame(α=Float64[], m=Int64[], converged=Bool[], free_energy=Float64[], fevals=Int64[])
    @showprogress for (a, m) in params
        updater = NGMRES(SD(a); m=m, warmup=0)
        scft = build_scft(updater, seedfile)
        ioconfig = my_ioconfig()
        scftconfig = SCFTConfig(max_iter=1000, tol=1e-8)
        config = JP.Config(scft=scftconfig, io=ioconfig)
        convergence = JP.solve!(scft, config)
        converged = (JP.Successful() == convergence)
        free_energy = JP.F(scft)
        push!(ngmres_sd_df, (a, m, converged, free_energy, scft.updater.evals[end]))
    end
    CSV.write(joinpath(csv_dir, "ngmres_sd_benchmark.csv"), ngmres_sd_df)
end
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

function sis_stat(csv_dir::String, seedfile::String)
    sis_df = DataFrame(α=Float64[], m=Int64[], converged=Bool[], free_energy=Float64[], fevals=Int64[])
    @showprogress for a in α
        updater = SIS(a)
        scft = build_scft(updater, seedfile)
        ioconfig = my_ioconfig()
        scftconfig = SCFTConfig(max_iter=3000, tol=1e-8)
        config = JP.Config(scft=scftconfig, io=ioconfig)
        convergence = JP.solve!(scft, config)
        converged = (JP.Successful() == convergence)
        free_energy = JP.F(scft) 
        push!(sis_df, (a, 0, converged, free_energy, scft.updater.evals[end]))    
    end
    CSV.write(joinpath(csv_dir, "sis_benchmark.csv"), sis_df)
end

function anderson_sis_stat(csv_dir::String, seedfile::String) 
    anderson_sis_df = DataFrame(α=Float64[], m=Int64[], converged=Bool[], free_energy=Float64[], fevals=Int64[])
    @showprogress for (a, m) in params
        updater = Anderson(SIS(a); m=m, warmup=0)
        scft = build_scft(updater, seedfile)
        ioconfig = my_ioconfig()
        scftconfig = SCFTConfig(max_iter=1000, tol=1e-8)
        config = JP.Config(scft=scftconfig, io=ioconfig)
        convergence = JP.solve!(scft, config)
        converged = (JP.Successful() == convergence)
        free_energy = JP.F(scft) 
        push!(anderson_sis_df, (a, m, converged, free_energy, scft.updater.evals[end]))    
    end
    CSV.write(joinpath(csv_dir, "anderson_sis_benchmark.csv"), anderson_sis_df)
end

function oaccel_sis_stat(csv_dir::String, seedfile::String)
    oaccel_sis_df = DataFrame(α=Float64[], m=Int64[], converged=Bool[], free_energy=Float64[], fevals=Int64[])
    @showprogress for (a, m) in params
        updater = OACCEL(SIS(a); m=m, warmup=0)
        scft = build_scft(updater, seedfile)
        ioconfig = my_ioconfig()
        scftconfig = SCFTConfig(max_iter=1000, tol=1e-8)
        config = JP.Config(scft=scftconfig, io=ioconfig)
        convergence = JP.solve!(scft, config)
        converged = (JP.Successful() == convergence)
        free_energy = JP.F(scft)
        push!(oaccel_sis_df, (a, m, converged, free_energy, scft.updater.evals[end]))
    end
    CSV.write(joinpath(csv_dir, "oaccel_sis_benchmark.csv"), oaccel_sis_df)
end

function ngmres_sis_stat(csv_dir::String, seedfile::String)
    ngmres_sis_df = DataFrame(α=Float64[], m=Int64[], converged=Bool[], free_energy=Float64[], fevals=Int64[])
    @showprogress for (a, m) in params
        updater = NGMRES(SIS(a); m=m, warmup=0)
        scft = build_scft(updater, seedfile)
        ioconfig = my_ioconfig()
        scftconfig = SCFTConfig(max_iter=1000, tol=1e-8)
        config = JP.Config(scft=scftconfig, io=ioconfig)
        convergence = JP.solve!(scft, config)
        converged = (JP.Successful() == convergence)
        free_energy = JP.F(scft)
        push!(ngmres_sis_df, (a, m, converged, free_energy, scft.updater.evals[end]))
    end
    CSV.write(joinpath(csv_dir, "ngmres_sis_benchmark.csv"), ngmres_sis_df)
end
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

function etd_stat(csv_dir::String, seedfile::String)
    etd_df = DataFrame(α=Float64[], m=Int64[], converged=Bool[], free_energy=Float64[], fevals=Int64[])
    @showprogress for a in α
        updater = ETD(a)
        scft = build_scft(updater, seedfile)
        ioconfig = my_ioconfig()
        scftconfig = SCFTConfig(max_iter=3000, tol=1e-8)
        config = JP.Config(scft=scftconfig, io=ioconfig)
        convergence = JP.solve!(scft, config)
        converged = (JP.Successful() == convergence)
        free_energy = JP.F(scft) 
        push!(etd_df, (a, 0, converged, free_energy, scft.updater.evals[end]))
    end
    CSV.write(joinpath(csv_dir, "etd_benchmark.csv"), etd_df)
end

function anderson_etd_stat(csv_dir::String, seedfile::String)
    anderson_etd_df = DataFrame(α=Float64[], m=Int64[], converged=Bool[], free_energy=Float64[], fevals=Int64[])
    @showprogress for (a, m) in params
        updater = Anderson(ETD(a); m=m, warmup=0)
        scft = build_scft(updater, seedfile)
        ioconfig = my_ioconfig()
        scftconfig = SCFTConfig(max_iter=1000, tol=1e-8)
        config = JP.Config(scft=scftconfig, io=ioconfig)
        convergence = JP.solve!(scft, config)
        converged = (JP.Successful() == convergence)
        free_energy = JP.F(scft)
        push!(anderson_etd_df, (a, m, converged, free_energy, scft.updater.evals[end]))
    end
    CSV.write(joinpath(csv_dir, "anderson_etd_benchmark.csv"), anderson_etd_df)
end

function oaccel_etd_stat(csv_dir::String, seedfile::String)
    oaccel_etd_df = DataFrame(α=Float64[], m=Int64[], converged=Bool[], free_energy=Float64[], fevals=Int64[])
    @showprogress for (a, m) in params
        updater = OACCEL(ETD(a); m=m, warmup=0)
        scft = build_scft(updater, seedfile)
        ioconfig = my_ioconfig()
        scftconfig = SCFTConfig(max_iter=1000, tol=1e-8)
        config = JP.Config(scft=scftconfig, io=ioconfig)
        convergence = JP.solve!(scft, config)
        converged = (JP.Successful() == convergence)
        free_energy = JP.F(scft)
        push!(oaccel_etd_df, (a, m, converged, free_energy, scft.updater.evals[end]))
    end
    CSV.write(joinpath(csv_dir, "oaccel_etd_benchmark.csv"), oaccel_etd_df)
end

function ngmres_etd_stat(csv_dir::String, seedfile::String)
    ngmres_etd_df = DataFrame(α=Float64[], m=Int64[], converged=Bool[], free_energy=Float64[], fevals=Int64[])
    @showprogress for (a, m) in params
        updater = NGMRES(ETD(a); m=m, warmup=0)
        scft = build_scft(updater, seedfile)
        ioconfig = my_ioconfig()
        scftconfig = SCFTConfig(max_iter=1000, tol=1e-8)
        config = JP.Config(scft=scftconfig, io=ioconfig)
        convergence = JP.solve!(scft, config)
        converged = (JP.Successful() == convergence)
        free_energy = JP.F(scft)
        push!(ngmres_etd_df, (a, m, converged, free_energy, scft.updater.evals[end]))
    end
    CSV.write(joinpath(csv_dir, "ngmres_etd_benchmark.csv"), ngmres_etd_df)
end

function one_seed_stat(csv_dir::String, seedfile::String)
    sd_stat(csv_dir, seedfile)
    anderson_sd_stat(csv_dir, seedfile)
    oaccel_sd_stat(csv_dir, seedfile)
    ngmres_sd_stat(csv_dir, seedfile)
    sis_stat(csv_dir, seedfile)
    anderson_sis_stat(csv_dir, seedfile)
    oaccel_sis_stat(csv_dir, seedfile)
    ngmres_sis_stat(csv_dir, seedfile)
    etd_stat(csv_dir, seedfile)
    anderson_etd_stat(csv_dir, seedfile)
    oaccel_etd_stat(csv_dir, seedfile)
    ngmres_etd_stat(csv_dir, seedfile)
end
