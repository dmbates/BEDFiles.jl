__precompile__()

module BEDFiles
    using Missings, Mmap
    import Base: IndexStyle, eltype, getindex, length, size

"""
    BEDvals

`Vector{Union{UInt8,Missings.Missing}}` of the possible values in a BEDFile
"""
    const bedvals = Union{UInt8,Missing}[0x00, missing, 0x01, 0x02]

"""
    BEDFile

Raw .bed file as a shared, memory-mapped Matrix{UInt8}.  The number of rows, `m`
is stored separately because it is not uniquely determined by the size of the `data` field.
"""
    struct BEDFile
        data::Matrix{UInt8}
        m::Int
    end
    function BEDFile(bnm::AbstractString) # bnm = base file name without extension
        m = countlines(string(bnm, ".fam"))
        n  = countlines(string(bnm, ".bim"))
        bednm = string(bnm, ".bed")
        bb = open(bednm, "r") do io
            read(io, UInt16) == 0x1b6c || throw(ArgumentError("wrong magic number in file $bednm"))
            isone(read(io, UInt8)) || throw(ArgumentError(".bed file, $bednm, is not in correct orientation"))
            Mmap.mmap(io, Matrix{UInt8}, ((m + 3) >> 2, n))
        end
        BEDFile(bb, m)
    end

"""
    BEDColumn

A single column from a `BEDFile`. Fields are `data`, a view of the underlying [`BEDFile`](@ref) column,
and `m`, the length of the column.
"""
    struct BEDColumn <: AbstractVector{Union{Missing,UInt8}}
        data::SubArray{UInt8,1,Array{UInt8,2},Tuple{Base.Slice{Base.OneTo{Int}},Int}}
        m::Int
    end    
    BEDColumn(f::BEDFile, j::Number) = BEDColumn(view(f.data, :, j), f.m)

    function Base.getindex(c::BEDColumn, i)
        0 < i ≤ c.m || throw(BoundsError("attempt to access $(c.m) element BEDColumn at index [$i]"))
        ip3 = i + 3
        @inbounds(bedvals[1 + ((c.data[ip3 >> 2] >> ((ip3 & 0x03) << 1)) & 0x03)])
    end

    Base.IndexStyle(::Type{BEDColumn}) = IndexLinear()

    function Base.iterate(c::BEDColumn, i=1)
        i ≤ c.m || return nothing
        ip3 = i + 3  # code repeated from getindex method b/c check on i can be skipped
        convert(Union{Missing,UInt8}, @inbounds(bedvals[1 + ((c.data[ip3 >> 2] >> ((ip3 & 0x03) << 1)) & 0x03)])), i + 1
    end    

    Base.length(c::BEDColumn) = c.m

    Base.eltype(::Type{BEDColumn}) = Union{Missing,UInt8}

"""
    outer!(pr::Matrix{Union{Missing,Int}}, c::BEDColumn)

Return `pr + c*c'`
"""
    function outer!(pr::Matrix{Union{Missing,Int}}, c::BEDColumn)
        j = 0
        for v in c
            prj = view(pr, :, j += 1)
            if ismissing(v)
                fill!(prj, missing)
            elseif isone(v)
                @. prj += c
            elseif !iszero(v)
                @. prj += v * c
            end
        end
        pr
    end

"""
    outer(c::BEDColumn)

Return `c * adjoint(c)`
"""
    outer(c::BEDColumn) = outer!(zeros(Union{Missing,Int}, (c.m, c.m)), c)

    Base.size(c::BEDColumn) = (c.m,)

    Base.size(f::BEDFile) = f.m, size(f.data, 2)
    Base.getindex(g::BEDFile, i::Integer, j::Integer) = BEDColumn(g, j)[i]

    export
        BEDFile,
        BEDColumn,

        bedvals,
        outer!,
        outer

end # module
