__precompile__()

module BEDFiles
    using LinearAlgebra, Mmap, OffsetArrays, SparseArrays, StaticArrays, Statistics, StatsBase
    export BEDFile, bedvals, counts, grm, maf, maf!, minorallele, minorallele!, 
        mean, mean!, missingpos, missingrate, missingrate!, std, var
    
    
    include("bedfile.jl")
    include("summarystats.jl")
    include("grm.jl")

"""
    BEDvals

`Vector{Union{UInt8,Missing}}` of the possible values in a BEDFile
"""
    const bedvals = OffsetArray{Union{Int8,Missing}}(undef, 0:3)
    bedvals[0] = 0
    bedvals[2] = 1
    bedvals[3] = 2

    datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end # module
