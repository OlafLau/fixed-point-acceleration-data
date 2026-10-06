using Pkg
if haskey(ENV, "POLYORDER_PROJECT")
    Pkg.activate(ENV["POLYORDER_PROJECT"])
end

using Polymer
using Scattering
using Polyorder
using Polyorder: SD, Anderson, SIS, ETD, NGMRES, OACCEL, BB
import Polyorder as JP
using DataFrames
using CSV
using Makie
using CairoMakie
using ProgressMeter
using Logging
using Statistics

function AB_lamellar_xN30(updater::SCFTAlgorithm)
    scft = let
        ds = 0.01
        uc = UnitCell(Square(), 4.4826)
        lat = BravaisLattice(uc)
        system = AB_system(χN=30.0, fA=0.50)
        NoncyclicChainSCFT(system, lat, ds, updater=updater)
    end
    return scft
end

function AB_cylinder_xN30(updater::SCFTAlgorithm)
    scft = let
        ds = 0.01
        uc = UnitCell(Square(), 3.951)
        lat = BravaisLattice(uc)
        system = AB_system(χN=30.0, fA=0.25)
        NoncyclicChainSCFT(system, lat, ds, updater=updater)
    end
    return scft
end

function chainAB3()
    
end

function AB3_system(fA::Float64=0.5, χN::Float64=30.0)
    sA = KuhnSegment(:A)
    sB = KuhnSegment(:B)
    eb = BranchPoint(:EB)
    A = PolymerBlock(:A, sA, fA, FreeEnd(:A), eb)
    B1 = PolymerBlock(:B1, sB, (1-fA)/3, eb, FreeEnd(:B1))
    B2 = PolymerBlock(:B2, sB, (1-fA)/3, eb, FreeEnd(:B2))
    B3 = PolymerBlock(:B3, sB, (1-fA)/3, eb, FreeEnd(:B2))
    block_copolymer = BlockCopolymer(:AB3, [A, B1, B2, B3])
    system = 
        let
            chain = block_copolymer  # Create a polymer chain
            polymer = Component(chain)  # Wrap the chain into a component
            # Create an interaction map between different species
            χNmap = Dict([:A, :B]=>χN,)
            PolymerSystem([polymer], χNmap)  # Create a polymer system
        end
    return system
end

# Gyroid phase
function AB3_gyr(system::AbstractSystem, updater::SCFTAlgorithm)
    scft = let
        ds = 0.01
        uc = UnitCell(Cubic(), 8.2813)
        lat = BravaisLattice(uc)
        NoncyclicChainSCFT(system, lat, ds, updater=updater)
    end
    return scft
end

function AB3_gyr(updater::SCFTAlgorithm)
    system = AB3_system(0.5, 20.0)
    scft = let
        ds = 0.01
        uc = UnitCell(Cubic(), 8.2813)
        lat = BravaisLattice(uc)
        NoncyclicChainSCFT(system, lat, ds, updater=updater)
    end
    return scft
end

# A15 phase
function AB3_sigma(system::AbstractSystem, updater::SCFTAlgorithm)
    scft = let
        ds = 0.01
        uc = UnitCell(Tetragonal(), 12.3042, 6.501)
        lat = BravaisLattice(uc)
        NoncyclicChainSCFT(system, lat, ds, updater=updater)
    end
    return scft
end

function AB3_sigma(updater::SCFTAlgorithm)
    system = AB3_system(0.34, 20.0)
    scft = let
        ds = 0.01
        uc = UnitCell(Tetragonal(), 12.3042, 6.501)
        lat = BravaisLattice(uc)
        NoncyclicChainSCFT(system, lat, ds, updater=updater)
    end
    return scft
end

# BCC phase
function AB3_bcc(system::AbstractSystem, updater::SCFTAlgorithm)
    scft = let
        ds = 0.01
        uc = UnitCell(Cubic(), 3.900)
        lat = BravaisLattice(uc)
        NoncyclicChainSCFT(system, lat, ds, updater=updater)
    end
    return scft
end

function AB3_bcc(updater::SCFTAlgorithm)
    system = AB3_system(0.4, 16.0)
    scft = let
        ds = 0.01
        uc = UnitCell(Cubic(), 3.900)
        lat = BravaisLattice(uc)
        NoncyclicChainSCFT(system, lat, ds, updater=updater)
    end

    return scft
end

# print all available functions
println(methods(AB_lamellar_xN30))
println(methods(AB_cylinder_xN30))
println(methods(AB3_system))
println(methods(AB3_gyr))
println(methods(AB3_sigma))
println(methods(AB3_bcc))




