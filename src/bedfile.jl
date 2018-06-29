"""
    BEDFile

Raw .bed file as a shared, memory-mapped Matrix{UInt8}.  The number of rows, `m`
is stored separately because it is not uniquely determined by the size of the `data` field.
"""
struct BEDFile
    data::Matrix{UInt8}
    m::Int
end
function BEDFile(bednm::AbstractString, m::Integer)
    data = open(bednm, "r") do io
        read(io, UInt16) == 0x1b6c || throw(ArgumentError("wrong magic number in file $bednm"))
        isone(read(io, UInt8)) || throw(ArgumentError(".bed file, $bednm, is not in correct orientation"))
        Mmap.mmap(io, Vector{UInt8})
    end
    drows = (m + 3) >> 2   # the number of rows in the Matrix{UInt8}
    n, r = divrem(length(data), drows)
    iszero(r) || throw(ArgumentError("filesize of $bednm is not a multiple of $drows"))
    BEDFile(reshape(data, (drows, n)), m)
end
BEDFile(nm::AbstractString) = BEDFile(nm, countlines(string(splitext(nm)[1], ".fam")))

"""
    columncounts(f::BEDFile)

Return a matrix of frequency counts for the columns
"""
function columncounts(f::BEDFile)
    m, n = size(f)
    counts = Matrix{Int}(undef, (4, n))
    Threads.@threads for j in 1:n
        counts!(view(counts, :, j), BEDColumn(f, j))
    end
    counts
end

Base.getindex(g::BEDFile, i::Integer, j::Integer) = BEDColumn(g, j)[i]

Base.size(f::BEDFile) = f.m, size(f.data, 2)

function Statistics.mean(f::BEDFile; dims)
    m, n = size(f)
    isone(dims) || throw(ArgumentError("mean(f::BEDFile; dims) only defined for dims = 1"))
    means = Matrix{Float64}(undef, (1, n))
    Threads.@threads for j in 1:n
        means[j] = mean(BEDColumn(f, j))
    end
    means
end
