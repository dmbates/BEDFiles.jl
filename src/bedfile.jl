"""
    BEDFile

Raw .bed file as a shared, memory-mapped Matrix{UInt8}.
- `data`: The compressed data matrix in which up to 4 SNP calls are stored in each `UInt8` element
- `columncounts`: a `4` by `n` `Int` matrix of column counts for each of the 4 possible calls
-
- `rowcounts`: a `4` by `m` array of row counts for each of the 4 possible calls
- `m`: the number of rows in the virtual array.

`m` is stored separately because it is not uniquely determined from the size of `data`.  If `data`
has `k` rows then `m` could be `4k - 3` or `4k - 2` or `4k - 1` or `4k` because the number of bytes
in a column of `data` is `k = (m + 3) ÷ 4`
"""
struct BEDFile <: AbstractMatrix{UInt8}
    data::Matrix{UInt8}
    columncounts::Matrix{Int}
    staticcounts::Base.ReinterpretArray{SArray{Tuple{4},Int,1,4},1,Int,Array{Int,1}}
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
    d, r = divrem(m, 4)
    ccounts = zeros(Int, (4, n))
    for j in 1:n
        offset = drows * (j - 1)
        for i in 1:d
            bb = data[offset + i]
            for s in 0:2:6
                ccounts[((bb >> s) & 0x03) + 1, j] += 1
            end
            if !iszero(r)
                bb = data[offset + d + 1]
                for s in 0:2:(2*(r-1))
                    ccounts[((bb >> s) & 0x03) + 1, j] += 1
                end
            end
        end
    end
    BEDFile(reshape(data, (drows, n)), ccounts, reinterpret(SVector{4, Int}, vec(ccounts)),
        zeros(Int, (4, m)), m)
end
BEDFile(nm::AbstractString, args...; kwargs...) = BEDFile(nm, countlines(string(splitext(nm)[1], ".fam")), args...; kwargs...)

function BEDFile(::UndefInitializer, m::Integer, n::Integer)
    BEDFile(Matrix{UInt8}(undef, ((m + 3) >> 2, n)), zeros(Int, (4, n)), zeros(Int, (4, m)), m)
end

function BEDFile(file::AbstractString, f::BEDFile)
    open(file, "w+") do io
        write(io, 0x1b6c)
        write(io, 0x01)
        write(io, f.data)
    end
    BEDFile(file, f.m, "r+")
end

function BEDFile(file::AbstractString, m::Integer, n::Integer)
    open(file, "w+") do io
        write(io, 0x1b6c)
        write(io, 0x01)
        write(io, fill(0x00, ((m + 3) >> 2, n)))
    end
    BEDFile(file, m, "r+")
end

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

function _copyto_additive!(v::AbstractVector{T}, f::BEDFile, j::Integer) where T <: AbstractFloat
    @inbounds for i in 1:f.m
        fij = f[i, j]
        v[i] = iszero(fij) ? zero(T) : isone(fij) ? T(NaN) : fij - 1
    end    
    v
end

function _copyto_dominant!(v::AbstractVector{T}, f::BEDFile, j::Integer) where T <: AbstractFloat
    @inbounds for i in 1:f.m
        fij = f[i, j]
        v[i] = iszero(fij) ? zero(T) : isone(fij) ? T(NaN) : 1
    end
    v
end

function _copyto_recessive!(v::AbstractVector{T}, f::BEDFile, j::Integer) where T <: AbstractFloat
    @inbounds for i in 1:f.m
        fij = f[i, j]
        v[i] = (iszero(fij) || fij == 2) ? zero(T) : isone(fij) ? T(NaN) : 1
    end    
    v
end

function Base.copyto!(
    v::AbstractVector{T}, 
    f::BEDFile, 
    j::Integer; 
    model::Symbol = :additive,
    center::Bool = false,
    scale::Bool = false,
    impute::Bool = false
    ) where T <: AbstractFloat
    if model == :additive
        _copyto_additive!(v, f, j)
    elseif model == :dominant
        _copyto_dominant!(v, f, j)
    elseif model == :recessive
        _copyto_recessive!(v, f, j)
    else
        throw(ArgumentError("model has to be :additive, :dominant, or :recessive; got $model"))
    end
    if center || scale || impute
        cc = _counts(f, 1)
        μ = model == :additive ? (cc[3, j] + 2cc[4, j]) / (cc[1, j] + cc[3, j] + cc[4, j]) : 
            model == :dominant ? (cc[3, j] +  cc[4, j]) / (cc[1, j] + cc[3, j] + cc[4, j]) :
            cc[4, j] / (cc[1, j] + cc[3, j] + cc[4, j])
        σ = model == :additive ? sqrt(μ * (1 - μ / 2)) : sqrt(μ * (1 - μ))
        doscale = scale && (σ > 0)
        @inbounds for i in 1:f.m
            impute && isnan(v[i]) && (v[i] = T(μ))
            center && (v[i] -= μ)
            doscale && (v[i] /= σ)
        end
    end
    v
end

function Base.copyto!(
    v::AbstractMatrix{T}, 
    f::BEDFile, 
    colinds::AbstractVector{<:Integer};
    kwargs...
    ) where T <: AbstractFloat
    for (vj, j) in enumerate(colinds)
        Base.copyto!(view(v, :, vj), f, j; kwargs...)
    end
    v
end

function Base.copyto!(
    v::AbstractMatrix{T}, 
    f::BEDFile, 
    colmask::AbstractVector{Bool};
    kwargs...
    ) where T <: AbstractFloat
    length(colmask) == size(f, 2) || throw(ArgumentError("`length(colmask)` does not match `size(f, 2)`"))
    vj = 1
    for j in 1:length(colmask)
        if colmask[j] 
            Base.copyto!(view(v, :, vj), f, j; kwargs...)
            vj += 1
        end
    end
    v
end

function Base.convert(t::Type{Vector{T}}, f::BEDFile, j::Integer; kwargs...) where T <: AbstractFloat
    Base.copyto!(Vector{T}(undef, f.m), f, j; kwargs...)
end
function Base.convert(t::Type{Matrix{T}}, f::BEDFile, colinds::AbstractVector{<:Integer}; kwargs...) where T <: AbstractFloat
    Base.copyto!(Matrix{T}(undef, f.m, length(colinds)), f, colinds; kwargs...)
end
function Base.convert(t::Type{Matrix{T}}, f::BEDFile, colmask::AbstractVector{Bool}; kwargs...) where T <: AbstractFloat
    Base.copyto!(Matrix{T}(undef, f.m, count(colmask)), f, colmask; kwargs...)
end
Base.convert(t::Type{Matrix{T}}, f::BEDFile; kwargs...) where T <: AbstractFloat = Base.convert(t, f, 1:size(f, 2); kwargs...)

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

"""
    grm(A; method=:GRM, maf_threshold=0.01)

Compute empirical kinship matrix from a BEDFile. Missing genotypes are imputed
on the fly by mean.

# Input  
- `f`: a BEDFile

# Optional Arguments
- `method`: `:GRM` (default), `:MoM`, or `Robust`
- `maf_threshold`: columns with MAF `<maf_threshold` are excluded; default 0.01
- `cinds`: indices or mask of columns to be used for calculating GRM
- `t`: Float type for calculating GRM
"""
function grm(
    f::BEDFile;
    method::Symbol = :GRM,
    maf_threshold::Real = 0.01,
    cinds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    t::Type{T} = Float64
    ) where T <: AbstractFloat
    mf = maf(f)
    colinds = something(cinds, mf .≥ maf_threshold)
    n = eltype(colinds) == Bool ? count(colinds) : length(colinds)
    G = Mmap.mmap(Matrix{t}, f.m, n)
    if method == :GRM
        Base.copyto!(G, f, colinds, model=:additive, impute=true, center=true, scale=true)
        Φ = G * transpose(G)
        Φ ./= 2n
    elseif method == :MoM
        Base.copyto!(G, f, colinds, model=:additive, impute=true)
        G .-= 1
        Φ = G * transpose(G)
        c = sum(x -> abs2(x) + abs2(1 - x), mf)
        shft, scal = n / 2 - c, 1 / (n - c)
        @inbounds @simd for i in eachindex(Φ)
            Φ[i] = (Φ[i] / 2 + shft) * scal
        end
    elseif method == :Robust
        Base.copyto!(G, f, colinds, model=:additive, center=true, impute=true)
        scal = sum(x -> 4x * (1 - x), mf)
        Φ = G * transpose(G)
        Φ ./= scal
    else
        throw(ArgumentError("method should be :GRM, :MoM, or :Robust; got $method"))
    end
    Φ
end # function grm

"""
    BEDFiles.filter(A[, min_success_rate_per_row, min_success_rate_per_col, maxiters])

Filter a BEDFile by genotyping success rate.

# Input
- `A`: a BEDFile.
- `min_success_rate_per_row`: threshold for SNP genotyping success rate.
- `min_success_rate_per_col`: threshold for person genotyping success rate.
- `maxiters`: maximum number of filtering iterations.

# Output
- `rmask`: BitVector indicating remaining rows.
- `cmask`: BitVector indicating remaining cols.
"""
function filter(
    f::BEDFile, 
    min_success_rate_per_row::Real = 0.98,
    min_success_rate_per_col::Real = 0.98,
    maxiters::Integer = 5)
    m, n = size(f)
    rc, cc = zeros(Int, m), zeros(Int, n)
    rmask, cmask = trues(m), trues(n)
    rmiss, cmiss = 1 - min_success_rate_per_row, 1 - min_success_rate_per_col
    for iter in 1:maxiters
        fill!(rc, 0)
        fill!(cc, 0)
        @inbounds for j in 1:n
            cmask[j] || continue
            for i in 1:m
                rmask[i] && f[i, j] == 0x01 && (rc[i] += 1; cc[j] += 1)
            end
        end
        rows, cols = count(rmask), count(cmask)
        @inbounds for j in 1:n
            cmask[j] = cmask[j] && cc[j] < cmiss * rows
        end
        @inbounds for i in 1:m
            rmask[i] = rmask[i] && rc[i] < rmiss * cols
        end
        count(cmask) == cols && count(rmask) == rows && break
        iter == maxiters && (@warn "success rate not satisfied; consider increase maxiters")
    end
    rmask, cmask
end

"""
    BEDFiles.filter(src, rowinds, colinds; des = src * ".filtered")

Filter `src` Plink files according to row indices `rowinds` and column indices 
`colinds` and write to a new set of Plink files `des`.

# Input
- `src`: source Plink file name without suffix ".bed", ".fam" or ".bim".
- `rowinds`: row indices.
- `colinds`: column indices.

# Keyword arguments
- `des`: output Plink file name; defualt is `src * ".filtered"`.
"""
function filter(
    src::AbstractString, 
    rowinds::AbstractVector{<:Integer},
    colinds::AbstractVector{<:Integer};
    des::AbstractString = src * ".filtered")
    # check source plink files
    isfile(src * ".bed") || throw(ArgumentError("$src.bed file not found"))
    isfile(src * ".bim") || throw(ArgumentError("$src.bim file not found"))
    isfile(src * ".fam") || throw(ArgumentError("$src.fam file not found"))
    # create row and column masks
    if eltype(rowinds) == Bool
        rmask = rowinds
    else
        rmask = falses(countlines(src * ".fam"))
        rmask[rowinds] .= true
    end
    if eltype(colinds) == Bool
        cmask = colinds
    else
        cmask = falses(countlines(src * ".bim"))
        cmask[colinds] .= true
    end
    m, n = count(rmask), count(cmask)
    bfsrc = BEDFile(src * ".bed")
    # write filtered bed file
    open(des * ".bed", "w+") do io
        write(io, 0x1b6c)
        write(io, 0x01)
        write(io, Matrix{UInt8}(undef, (m + 3) >> 2, n))
    end
    bfdes = BEDFile(des * ".bed", m, "r+")
    bfdes .= @view bfsrc[rmask, cmask]
    # write filtered fam file
    open(des * ".fam", "w") do io
        for (i, line) in enumerate(eachline(src * ".fam"))
            rmask[i] && println(io, line)
        end
    end
    # write filtered bim file
    open(des * ".bim", "w") do io
        for (j, line) in enumerate(eachline(src * ".bim"))
            cmask[j] && println(io, line)
        end
    end
    # output BEDFile
    bfdes
end
