StatsBase.counts(f::BEDFile; dims=:) = _counts(f, dims)

function _counts(f::BEDFile, dims::Integer)
    if dims == 1
        f.columncounts
    elseif dims == 2
        f.rowcounts
    else
        throw(ArgumentError("counts(f::BEDFile, dims=k) only defined for k = 1 or 2"))
    end
end

_counts(f::BEDFile, ::Colon) = sum(f.sccounts)

@inline nnonmiss(v::SVector{4,Int}) = v[1] + v[3] + v[4]

meanfromcounts(v::SVector{4,Int}) = (v[3] + 2v[4]) / nnonmiss(v)

Statistics.mean(f::BEDFile; dims=:) = _mean(f, dims)

function _mean(f::BEDFile,  dims::Integer)
    if dims == 1
        reshape(meanfromcounts.(f.sccounts), (1, size(f, 2)))
    elseif dims == 2
        reshape(meanfromcounts.(f.srcounts), (size(f, 1), 1))
    else
        throw(ArgumentError("mean(f::BEDFile, dims=k) only defined for k = 1 or 2"))
    end
end

_mean(f::BEDFile, ::Colon) = meanfromcounts(sum(f.sccounts))

function Statistics.mean!(out::AbstractMatrix{<:AbstractFloat}, f::BEDFile)
    k, l = size(out)
    m, n = size(f)
    errormsg = "size(out) should be (1,1) or (1,size(f,2)) or (size(f,1),1)"
    if k == 1
        if l == 1
            out[1] = _mean(f, :)
        elseif l == n
            map!(meanfromcounts, out, f.sccounts)
        else
            throw(ArgumentError(errormsg))
        end
    elseif k == m && l == 1
        map!(meanfromcounts, out, f.srcounts)
    else
        throw(ArgumentError(errmsg))
    end
    out
end

"""
    missingrate!(out::AbstractMatrix{<:AbstractFloat}, f::BEDFile)

Overwrite the contents of `out` with the missing rates according to the size of `out`.

If `size(out) == (1,size(f,2))` it is overwritten with the column-wise missing rates.
If `size(out) == (size(f,1),1)` it is overwritten with the row-wise missing rates.
If `size(out) == (1,1)` it is overwritten with the overall rate of missing values.
"""
function missingrate!(out::AbstractMatrix{<:AbstractFloat}, f::BEDFile)
    k, l = size(out)
    m, n = size(f)
    errormsg = "size(out) should be (1,1) or (1,size(f,2)) or (size(f,1),1)"
    if k == 1
        if l == 1
            out[1] = sum(v -> v[2], f.sccounts) / length(f)
        elseif l == n
            map!(v -> v[2] / m, out, f.sccounts)
        else
            throw(ArgumentError(errormsg))
        end
    elseif k == m && l == 1
        map!(v -> v[2] / n, out, f.srcounts)
    else
        throw(ArgumentError(errormsg))
    end
    out
end

"""
    missingrate(f::BEDFile; dims=:)

Return the rate of missing values in `f` by column (`dims=1`), by row (`dims=2`) or overall (`dims=:`, the default).
"""
missingrate(f::BEDFile; dims=:) = _missingrate(f, dims)

function _missingrate(f::BEDFile, dims::Integer)
    m, n = size(f)
    if dims == 1
        reshape([v[2] / m for v in f.sccounts], (1, n))
    elseif dims == 2
        reshape([v[2] / n for v in f.srcounts], (m, 1))
    else
        throw(ArgumentError("missingrate(f::BEDFile, dims=k) only defined for k = 1 or 2"))
    end
end

_missingrate(f::BEDFile, ::Colon) = sum(view(f.rowcounts, 2, :)) / length(f)

Statistics.var(f::BEDFile; corrected::Bool=true, mean=nothing, dims=:) = 
    _var(f, corrected, mean, dims)

function varfromcounts(v::SVector{4,Int}, corrected::Bool, mn)
    nnmiss = nnonmiss(v)
    if mn == nothing
        mn = (v[3] + 2v[4]) / nnmiss
    end
    (abs2(mn)*v[1] + abs2(1 - mn)*v[3] + abs2(2 - mn)*v[4]) / (nnmiss - Int(corrected))
end    

function _var(f::BEDFile, corrected::Bool, mean, dims::Integer)
    m, n = size(f)
    if dims == 1
        reshape([varfromcounts(v, corrected, mean) for v in f.sccounts], (1, n))
    elseif dims == 2
        reshape([varfromcounts(v, corrected, mean) for v in f.srcounts], (m, 1))
    else
        throw(ArgumentError("var(f::BEDFile, dims=k) only defined for k = 1 or 2"))
    end
end

function maffromcounts(v::SVector{4, Int})
    freq = (v[3] + 2v[4]) / (2 * nnonmiss(v))
    freq ≤ 0.5 ? freq : 1 - freq
end

maf!(out::AbstractVector{<:AbstractFloat}, f::BEDFile) = map!(maffromcounts, out, f.sccounts)

"""
    maf(f::BEDFile)

Return a vector of minor allele frequencies for the columns of `f`.

By definition the minor allele frequency is between 0 and 0.5
"""
maf(f::BEDFile) = maffromcounts.(f.sccounts)

minorallelefromcounts(v::SVector{4, Int}) = v[1] > v[4]

minorallele!(out::AbstractVector{Bool}, f::BEDFile) =
    map!(minorallelefromcounts, out, f.sccounts)

"""
    minorallele(f::BEDFile)

Return a `Vector{Bool}` indicating if the minor allele in each column is A2
"""
minorallele(f::BEDFile) = minorallelefromcounts.(f.sccounts)
