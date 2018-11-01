StatsBase.counts(f::BEDFile; dims=:) = _counts(f, dims)

function _counts(f::BEDFile, dims::Integer)
    isone(dims) && return f.columncounts
    if dims == 2
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

_counts(f::BEDFile, ::Colon) = sum(f.staticcounts)

Statistics.mean(f::BEDFile; dims=:) = _mean(f, dims)

function _mean(f::BEDFile,  dims::Integer)
    if isone(dims)
        mean!(Matrix{Float64}(undef, (1, size(f, 2))), f)
    elseif dims == 2
        mean!(Matrix{Float64}(undef, (size(f, 1), 1)), f)
    else
        throw(ArgumentError("mean(f::BEDFile, dims=k) only defined for k = 1 or 2"))
    end
end

nnonmiss(v::SVector{4,Int}) = v[1] + v[3] + v[4]

meanfromcounts(v::SVector{4,Int}) = (v[3] + 2v[4]) / nnonmiss(v)

_mean(f::BEDFile, ::Colon) = meanfromcounts(sum(f.staticcounts))

function Statistics.mean!(out::AbstractMatrix{<:AbstractFloat}, f::BEDFile)
    k, l = size(out)
    m, n = size(f)
    errormsg = "size(out) = $(size(out)) should be (1,1) or (1,$n) or ($m,1)"
    if isone(k)
        if isone(l)
            out[1] = _mean(f, :)
        elseif l == n
            map!(meanfromcounts, out, f.staticcounts)
        else
            throw(ArgumentError(errormsg))
        end
    elseif k == m && isone(l)
        rc = _counts(f, 2)
        @inbounds for i in 1:m
            out[i, 1] = (rc[3, i] + 2rc[4, i]) / (rc[1, i] + rc[3, i] + rc[4, i])
        end
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
    errormsg = "size(out) = $(size(out)) should be (1,1) or (1,$n) or ($m,1)"
    if isone(k)
        if isone(l)
            out[1] = sum(v -> v[2], f.staticcounts) / length(f)
        elseif l == n
            map!(v -> v[2] / m, out, f.staticcounts)
        else
            throw(ArgumentError(errormsg))
        end
    elseif k == m && isone(l)
        rc = _counts(f, 2)
        @inbounds for i in 1:m
            out[i] = rc[2, i] / n
        end
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
    if isone(dims)
        missingrate!(Matrix{Float64}(undef, (1, size(f, 2))), f)
    elseif dims == 2 
        missingrate!(Matrix{Float64}(undef, (f.m, 1)), f)
    else
        throw(ArgumentError("missingrate(f::BEDFile, dims=k) only defined for k = 1 or 2"))
    end
end

_missingrate(f::BEDFile, ::Colon) = sum(v -> v[2], f.staticcounts) / length(f)

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
    if isone(dims)
        reshape([varfromcounts(v, corrected, mean) for v in f.staticcounts], (1, size(f,2)))
    elseif dims == 2
        means = something(mean, Statistics.mean(f, dims=2))
        rc = _counts(f, 2)
        vars = Matrix{Float64}(undef, (m, 1))
        for i in 1:m
            mni = means[i]
            vars[i] = (abs2(mni)*rc[1,i] + abs2(1.0 - mni)*rc[3,i] + abs2(2.0 - mni)*rc[4,i]) /
                (rc[1,i] + rc[3,i] + rc[4,i] - Int(corrected))
        end
        vars
    else
        throw(ArgumentError("var(f::BEDFile, dims=k) only defined for k = 1 or 2"))
    end
end

function maffromcounts(v::SVector{4, Int})
    freq = (v[3] + 2v[4]) / (2 * nnonmiss(v))
    freq ≤ 0.5 ? freq : 1 - freq
end

maf!(out::AbstractVector{<:AbstractFloat}, f::BEDFile) = map!(maffromcounts, out, f.staticcounts)

"""
    maf(f::BEDFile)

Return a vector of minor allele frequencies for the columns of `f`.

By definition the minor allele frequency is between 0 and 0.5
"""
maf(f::BEDFile) = maf!(Vector{Float64}(undef, size(f, 2)), f)

minorallelefromcounts(v::SVector{4, Int}) = v[1] > v[4]

minorallele!(out::AbstractVector{Bool}, f::BEDFile) =
    map!(minorallelefromcounts, out, f.staticcounts)

"""
    minorallele(f::BEDFile)

Return a `Vector{Bool}` indicating if the minor allele in each column is A2
"""
minorallele(f::BEDFile) = minorallele!(Vector{Bool}(undef, size(f, 2)), f)
