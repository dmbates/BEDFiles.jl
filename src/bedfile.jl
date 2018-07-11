"""
    BEDFile

Raw .bed file as a shared, memory-mapped Matrix{UInt8}.  The number of rows, `m`
is stored separately because it is not uniquely determined by the size of the `data` field.
"""
struct BEDFile <: AbstractMatrix{UInt8}
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
columncounts(f::BEDFile) = colcounts!(zeros(Int, (4, size(f, 2))), f)

function colcounts!(counts, f::BEDFile)
    m, n = size(f)
    @inbounds for j in 1:n
        for i in 1:m
            counts[f[i, j] + 1, j] += 1
        end
    end    
    counts
end    

function Base.getindex(f::BEDFile, i::Int)  # Linear indexing
    d, r = divrem(i, f.m)
    f[r + 1, d + 1]
end
function Base.getindex(f::BEDFile, i::Integer, j::Integer)
    ip3 = i + 3
    (f.data[ip3 >> 2, j] >> ((ip3 & 0x03) << 1)) & 0x03
end

Base.eltype(f::BEDFile) = UInt8

Base.length(f::BEDFile) = f.m * size(f.data, 2)

Base.size(f::BEDFile) = f.m, size(f.data, 2)

Base.size(f::BEDFile, k::Integer) = 
    k == 1 ? f.m : k == 2 ? size(f.data, 2) : k > 2 ? 1 : error("Dimension out of range")

"""
    outer(f::BEDFile, colinds)
    outer(f::BEDFile)

Return the "outer product", `f * f'` using the `Float32[0, NaN, 1, 2]` encoding of `f`

The `colinds` argument, when given, causes the operation to be performed on that subset
of the columns.
"""
function outer(f::BEDFile, colinds::AbstractVector{<:Integer})
    m = size(f, 1)
    outer!(Symmetric(zeros(Float32, (m, m))), f, colinds)
end
outer(f::BEDFile) = outer(f, 1:size(f, 2))

function Base.copyto!(v::AbstractVector{T}, f::BEDFile, j) where T <: AbstractFloat
    for i in 1:f.m
        fij = f[i, j]
        v[i] = iszero(fij) ? zero(T) : isone(fij) ? T(NaN) : fij - 1
    end
    v
end

"""
    outer!(sy::Symmetric, f::BEDFile, colinds)

update `sy` with the sum of the outer products of the columns in `colind` from `f`
"""
function outer!(sy::Symmetric{T}, f::BEDFile, colinds::AbstractVector{<:Integer}) where T
    tempv = Vector{T}(undef, f.m)
    for j in colinds
        LinearAlgebra.BLAS.syr!(sy.uplo, one(T), copyto!(tempv, f, j), sy.data)
    end
    sy
end

function Statistics.mean(f::BEDFile; dims)
    m, n = size(f)
    isone(dims) || throw(ArgumentError("mean(f::BEDFile; dims) only defined for dims = 1"))
    means = Matrix{Float64}(undef, (1, n))
    @inbounds for j in 1:n
        s = 0.0
        nnmiss = 0
        for i in 1:m
            fij = f[i, j]
            if fij ≠ 1
                nnmiss += 1
            end
            if fij > 1
                s += fij - 1
            end
        end
        means[j] = s/nnmiss
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
