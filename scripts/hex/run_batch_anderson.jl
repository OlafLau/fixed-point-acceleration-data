include(joinpath(@__DIR__, "base.jl"))

anderson_sd_stat = Dict()
warmup_list = [w for w in 0:10:200]
for warmup in warmup_list
    updater = Anderson(SD(0.4); m=30, warmup=warmup)
    one_stat_result = run_batch_stat(updater)
    anderson_sd_stat[string(warmup)] = one_stat_result
end
# save stat_result to JLD2 file
outfile = joinpath(STAT_DIR, "anderson_sd0.4_stat.jld2")
@save outfile anderson_sd_stat

# --------------------------------------------------------------------------------------------------
anderson_sis_stat = Dict()
warmup_list = [w for w in 0:10:200]
for warmup in warmup_list
    updater = Anderson(SIS(0.4); m=30, warmup=warmup)
    one_stat_result = run_batch_stat(updater)
    anderson_sis_stat[string(warmup)] = one_stat_result
end
# save stat_result to JLD2 file
outfile = joinpath(STAT_DIR, "anderson_sis0.4_stat.jld2")
@save outfile anderson_sis_stat

# --------------------------------------------------------------------------------------------------
anderson_etd_stat = Dict()
warmup_list = [w for w in 0:10:200]
for warmup in warmup_list
    updater = Anderson(ETD(0.4); m=30, warmup=warmup)
    one_stat_result = run_batch_stat(updater)
    anderson_etd_stat[string(warmup)] = one_stat_result
end
outfile = joinpath(STAT_DIR, "anderson_etd0.4_stat.jld2")
@save outfile anderson_etd_stat
