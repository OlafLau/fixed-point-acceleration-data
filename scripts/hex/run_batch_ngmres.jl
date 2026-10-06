include(joinpath(@__DIR__, "base.jl"))

warmup_list = [w for w in 0:10:200]
ngmres_sd_stat = Dict()

for warmup in warmup_list
    updater = NGMRES(SD(0.4); m=30, warmup=warmup)
    one_stat_result = run_batch_stat(updater)
    ngmres_sd_stat[string(warmup)] = one_stat_result
end
# save stat_result to JLD2 file
outfile = joinpath(STAT_DIR, "ngmres_sd0.4_stat.jld2")
@save outfile ngmres_sd_stat

# --------------------------------------------------------------------------------------------------
ngmres_sis_stat = Dict()
for warmup in warmup_list
    updater = NGMRES(SIS(0.4); m=30, warmup=warmup)
    one_stat_result = run_batch_stat(updater)
    ngmres_sis_stat[string(warmup)] = one_stat_result
end
# save stat_result to JLD2 file
outfile = joinpath(STAT_DIR, "ngmres_sis0.4_stat.jld2")
@save outfile ngmres_sis_stat

# --------------------------------------------------------------------------------------------------
ngmres_etd_stat = Dict()
for warmup in warmup_list
    updater = NGMRES(ETD(0.4); m=30, warmup=warmup)
    one_stat_result = run_batch_stat(updater)
    ngmres_etd_stat[string(warmup)] = one_stat_result
end
outfile = joinpath(STAT_DIR, "ngmres_etd0.4_stat.jld2")
@save outfile ngmres_etd_stat
