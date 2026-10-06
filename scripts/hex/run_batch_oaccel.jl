include(joinpath(@__DIR__, "base.jl"))

oaccel_sd_stat = Dict()
warmup_list = [w for w in 0:10:200]
for warmup in warmup_list
    updater = OACCEL(SD(0.4); m=30, warmup=warmup)
    one_stat_result = run_batch_stat(updater)
    oaccel_sd_stat[string(warmup)] = one_stat_result
end
# save stat_result to JLD2 file
outfile = joinpath(STAT_DIR, "oaccel_sd0.4_stat.jld2")
@save outfile oaccel_sd_stat

# --------------------------------------------------------------------------------------------------
oaccel_sis_stat = Dict()
warmup_list = [w for w in 0:10:200]
for warmup in warmup_list
    updater = OACCEL(SIS(0.4); m=30, warmup=warmup)
    one_stat_result = run_batch_stat(updater)
    oaccel_sis_stat[string(warmup)] = one_stat_result
end
# save stat_result to JLD2 file
outfile = joinpath(STAT_DIR, "oaccel_sis0.4_stat.jld2")
@save outfile oaccel_sis_stat

# --------------------------------------------------------------------------------------------------
oaccel_etd_stat = Dict()
warmup_list = [w for w in 0:10:200]
for warmup in warmup_list
    updater = OACCEL(ETD(0.4); m=30, warmup=warmup)
    one_stat_result = run_batch_stat(updater)
    oaccel_etd_stat[string(warmup)] = one_stat_result
end
outfile = joinpath(STAT_DIR, "oaccel_etd0.4_stat.jld2")
@save outfile oaccel_etd_stat
