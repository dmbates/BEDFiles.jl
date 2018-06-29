
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
    @inbounds(((c.data[ip3 >> 2] >> ((ip3 & 0x03) << 1)) & 0x03))
end

Base.IndexStyle(::Type{BEDColumn}) = IndexLinear()

StatsBase.counts(c::BEDColumn) = counts!(Vector{Int}(undef, 4), c)

function Base.iterate(c::BEDColumn, (i, b)=(1, 0x00))
    i ≤ c.m || return nothing
    ip3 = i + 3
    s = (ip3 & 0x03) << 1 # == 2*rem(ip3, 4) - the size of the bitshift
    iszero(s) && @inbounds(b = c.data[ip3 >> 2])
    (b >> s) & 0x03, (i + 1, b)
end    

Base.length(c::BEDColumn) = c.m

Base.eltype(::Type{BEDColumn}) = UInt8

Base.size(c::BEDColumn) = (c.m,)

"""
    sumnnmiss(c::BEDColumn)

Return the sum of the column in the [0, missing, 1, 2] encoding and the number of nonmissing values
"""
function sumnnmiss(c::BEDColumn)
    nnmiss = 0 # number of non-missing
    s = 0
    for v in c
        if v ≠ 1
            nnmiss += 1
        end
        if v > 1
            s += v - 1
        end
    end
    s, nnmiss
end

function Statistics.mean(c::BEDColumn)
    s, nnmiss = sumnnmiss(c)
    s / nnmiss
end


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
    counts!(counts::Vector, c::BEDColumn)

Return `counts` overwritten with the counts of the raw values in `c`
"""
function counts!(counts::AbstractVector{<:Integer}, c::BEDColumn)
    length(counts) ≥ 4 || throw(ArgumentError("length(counts) = $(length(counts)) should be at least 4"))
    fill!(counts, 0)
    for v in c
        @inbounds counts[v + 1] += 1
    end
    counts
end


