__precompile__()

module BEDFiles
    using LinearAlgebra, Missings, Mmap, OffsetArrays, SparseArrays, Statistics, StatsBase
    import Base: IndexStyle, convert, copyto!, eltype, getindex, length, size
    import Statistics: mean
    import StatsBase: counts
    export BEDFile, BEDColumn, bedvals, columncounts, counts, mean, missingpos, outer, outer!
    
    include("bedfile.jl")
    include("bedcolumn.jl")

"""
    BEDvals

`Vector{Union{UInt8,Missings.Missing}}` of the possible values in a BEDFile
"""
    const bedvals = OffsetArray{Union{Int8,Missing}}(undef, 0:3)
    bedvals[0] = 0
    bedvals[2] = 1
    bedvals[3] = 2

end # module
