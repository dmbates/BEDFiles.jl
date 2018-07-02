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
        read(io, UInt8) == 0x01 || throw(ArgumentError(".bed file, $bednm, is not in correct orientation"))
        Mmap.mmap(io)
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

"""
    outer(f::BEDFile, colinds)
    outer(f::BEDFile)

Return the "outer product", `f * f'` using the `Float32[0, NaN, 1, 2]` encoding of `f`

The `colinds` argument, when given, causes the operation to be performed on that subset
of the columns.
"""
function outer(f::BEDFile, colinds::AbstractVector{<:Integer})
    m, n = size(f)
    outer!(Symmetric(zeros(Float32, (m,m))), f, colinds)
end
outer(f::BEDFile) = outer(f, 1:size(f)[2])

"""
    outer!(sy::Symmetric, f::BEDFile, colinds)

update `sy` with the sum of the outer products of the columns in `colind` from `f`
"""
function outer!(sy::Symmetric{T}, f::BEDFile, colinds::AbstractVector{<:Integer}) where T <: AbstractFloat
    m, n = size(f)
    tempv = Vector{T}(undef, m)
    for j in colinds
        LinearAlgebra.BLAS.syr!(sy.uplo, 1.0f0, copyto!(tempv, BEDColumn(f, j)), sy.data)
    end
    sy
end

Base.size(f::BEDFile) = f.m, size(f.data, 2)

function Statistics.mean(f::BEDFile; dims)
    m, n = size(f)
    isone(dims) || throw(ArgumentError("mean(f::BEDFile; dims) only defined for dims = 1"))
    means = Matrix{Float64}(undef, (1, n))
    Threads.@threads for j in 1:n
        @inbounds means[j] = mean(BEDColumn(f, j))
    end
    means
end

"""
    missingpos(f::BEDFile)

Return a `SparseMatrixCSC` of the same size as `f` indicating the positions with missing data
"""
function missingpos(f::BEDFile)
    m, n = size(f)
    colptr = sizehint!(Int32[1], n)
    rowval = Int32[]
    for j in 1:n
        msngpos = findall(isone.(BEDColumn(f, j)))
        append!(rowval, msngpos)
        push!(colptr, colptr[end] + length(msngpos))
    end
    SparseMatrixCSC(m, n, colptr, rowval, ones(Int8, length(rowval)))
end
