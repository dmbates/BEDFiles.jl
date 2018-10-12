"""
    BEDFile

Raw .bed file as a shared, memory-mapped Matrix{UInt8}.  The number of rows, `m`
is stored separately because it is not uniquely determined by the size of the `data` field.
"""
struct BEDFile <: AbstractMatrix{UInt8}
    data::Matrix{UInt8}
    columncounts::Matrix{Int}
    rowcounts::Matrix{Int}
    m::Int
end
function BEDFile(bednm::AbstractString, m::Integer, args...; kwargs...)
    data = open(bednm, args...; kwargs...) do io
        read(io, UInt16) == 0x1b6c || throw(ArgumentError("wrong magic number in file $bednm"))
        read(io, UInt8) == 0x01 || throw(ArgumentError(".bed file, $bednm, is not in correct orientation"))
        Mmap.mmap(io)
    end
    drows = (m + 3) >> 2   # the number of rows in the Matrix{UInt8}
    n, r = divrem(length(data), drows)
    iszero(r) || throw(ArgumentError("filesize of $bednm is not a multiple of $drows"))
    BEDFile(reshape(data, (drows, n)), zeros(Int, (4, n)), zeros(Int, (4, m)), m)
end
BEDFile(nm::AbstractString, args...; kwargs...) = BEDFile(nm, countlines(string(splitext(nm)[1], ".fam")), args...; kwargs...)

StatsBase.counts(f::BEDFile; dims=:) = _counts(f, dims)

function _counts(f::BEDFile, dims::Integer)
    if isone(dims)
        cc = f.columncounts
        if all(iszero, cc)
            m, n = size(f)
            @inbounds for j in 1:n
                for i in 1:m
                    cc[f[i, j] + 1, j] += 1            
                end
            end
        end
        return cc
    elseif dims == 2
        rc = f.rowcounts
        if all(iszero, rc)
            m, n = size(f)
            @inbounds for j in 1:n
                for i in 1:m
                    rc[f[i, j] + 1, i] += 1            
                end
            end
        end
        return rc
    else
        throw(ArgumentError("counts(f::BEDFile, dims=k) only defined for k = 1 or 2"))
    end
end

_counts(f::BEDFile, ::Colon) = sum(_counts(f, 1), dims=2)

function Base.getindex(f::BEDFile, i::Int)  # Linear indexing
    d, r = divrem(i - 1, f.m)
    f[r + 1, d + 1]
end

@inline function Base.getindex(f::BEDFile, i::Integer, j::Integer)
    @boundscheck checkbounds(f, i, j)
    ip3 = i + 3
    (f.data[ip3 >> 2, j] >> ((ip3 & 0x03) << 1)) & 0x03
end

function Base.setindex!(f::BEDFile, x::UInt8, i::Int)  # Linear indexing
    d, r = divrem(i - 1, f.m)
    Base.setindex!(f, x, r + 1, d + 1)
end

@inline function Base.setindex!(f::BEDFile, x::UInt8, i::Integer, j::Integer)
    @boundscheck checkbounds(f, i, j)
    ip3 = i + 3
    shft = (ip3 & 0x03) << 1
    mask = ~(0x03 << shft)
    f.data[ip3 >> 2, j] = (f.data[ip3 >> 2, j] & mask) | (x << shft)
    x
end

Base.eltype(f::BEDFile) = UInt8

Base.length(f::BEDFile) = f.m * size(f.data, 2)

Statistics.mean(f::BEDFile; dims=:) = _mean(f, dims)

function _mean(f::BEDFile,  dims::Integer)
    m, n = size(f)
    if isone(dims)
        cc = _counts(f, 1)   # need to use extractor to force evaluation if needed
        means = Matrix{Float64}(undef, (1, n))
        @inbounds for j in 1:n
            means[j] = (cc[3, j] + 2*cc[4, j]) / (cc[1, j] + cc[3, j] + cc[4, j])
        end
        return means
    elseif dims == 2
        rc = _counts(f, 2)
        means = Matrix{Float64}(undef, (m, 1))
        @inbounds for i in 1:m
            means[i] = (rc[3, i] + 2*rc[4, i]) / (rc[1, i] + rc[3, i] + rc[4, i])
        end
        return means
    else
        throw(ArgumentError("mean(f::BEDFile, dims=k) only defined for k = 1 or 2"))
    end
end

function _mean(f::BEDFile, ::Colon)
    rc = _counts(f, 2)
    (sum(view(rc, 3, :)) + 2*sum(view(rc, 4, :))) / sum(view(rc, [1, 3, 4], :))
end

Base.size(f::BEDFile) = f.m, size(f.data, 2)

Base.size(f::BEDFile, k::Integer) = 
    k == 1 ? f.m : k == 2 ? size(f.data, 2) : k > 2 ? 1 : error("Dimension out of range")

Statistics.var(f::BEDFile; corrected::Bool=true, mean=nothing, dims=:) = _var(f, corrected, mean, dims)

function _var(f::BEDFile, corrected::Bool, mean, dims::Integer)
    m, n = size(f)
    means = something(mean, Statistics.mean(f, dims=dims))
    if isone(dims)
        cc = _counts(f, 1)
        vars = Matrix{Float64}(undef, (1, n))
        for j in 1:n
            mnj = means[j]
            vars[j] = (abs2(mnj)*cc[1,j] + abs2(1.0 - mnj)*cc[3,j] + abs2(2.0 - mnj)*cc[4,j]) /
                (cc[1,j] + cc[3,j] + cc[4,j] - (corrected ? 1 : 0))
        end
        return vars
    elseif dims == 2
        rc = _counts(f, 2)
        vars = Matrix{Float64}(undef, (m, 1))
        for i in 1:m
            mni = means[i]
            vars[i] = (abs2(mni)*rc[1,i] + abs2(1.0 - mni)*rc[3,i] + abs2(2.0 - mni)*rc[4,i]) /
                (rc[1,i] + rc[3,i] + rc[4,i] - (corrected ? 1 : 0))
        end
        return vars
    end
    throw(ArgumentError("var(f::BEDFile, dims=k) only defined for k = 1 or 2"))
end

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

"""
    missingpos(f::BEDFile)

Return a `SparseMatrixCSC{Bool,Int32}` of the same size as `f` indicating the positions with missing data
"""
function missingpos(f::BEDFile)
    m, n = size(f)
    colptr = sizehint!(Int32[1], n + 1)
    rowval = Int32[]
    @inbounds for j in 1:n
        msngpos = Int32[]
        for i in 1:m
            isone(f[i, j]) && push!(msngpos, i)
        end
        append!(rowval, msngpos)
        push!(colptr, colptr[end] + length(msngpos))
    end
    SparseMatrixCSC(m, n, colptr, rowval, fill(true, length(rowval)))
end
