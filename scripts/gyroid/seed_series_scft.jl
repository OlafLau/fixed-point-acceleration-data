include(joinpath(@__DIR__, "base.jl"))
include(joinpath(@__DIR__, "stat_base.jl"))

# xN = [20.0]
# fA = [0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.54, 0.55]

# xN = [20.0]
# fA = [0.51, 0.53]

xN = [20.0]
fA = [0.52]
xN_fA = [(chi, f) for chi in xN for f in fA]

for (xN, fA) in xN_fA
    seedfile = joinpath(@__DIR__, "seeds", "xN$(xN)_fA$(fA)_fields.h5")
    csv_dir = joinpath(@__DIR__, "results", "xN$(xN)_fA$(fA)_csv")
    mkpath(csv_dir)
    one_seed_stat(csv_dir, seedfile)
end
