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

_counts(f::BEDFile, ::Colon) = sum(f.staticcounts, dims=2)

Statistics.mean(f::BEDFile; dims=:) = _mean(f, dims)

meanfromcounts(v::SVector{4,Int}) = (v[3] + 2v[4]) / (v[1] + v[3] + v[4])

function _mean(f::BEDFile,  dims::Integer)
    if isone(dims)
        reshape(meanfromcounts.(f.staticcounts), (1, size(f, 2)))
    elseif dims == 2
        rc = _counts(f, 2)
        m = size(f, 1)
        means = Matrix{Float64}(undef, (m, 1))
        @inbounds for i in 1:m
            means[i] = (rc[3, i] + 2*rc[4, i]) / (rc[1, i] + rc[3, i] + rc[4, i])
        end
        means
    else
        throw(ArgumentError("mean(f::BEDFile, dims=k) only defined for k = 1 or 2"))
    end
end

function _mean(f::BEDFile, ::Colon)
    v = sum(f.staticcounts)
    (v[3] + 2v[4]) / (v[1] + v[3] + v[4])
end

function Statistics.mean!(out::AbstractMatrix{<:AbstractFloat}, f::BEDFile)
    k, l = size(out)
    m, n = size(f)
    if k == 1 && l == n
        map!(meanfromcounts, out, f.staticcounts)
    elseif k == m && l == 1
        rc = _counts(f, 2)
        @inbounds for i in 1:m
            out[i, 1] = (rc[3, i] + 2rc[4, i]) / (rc[1, i] + rc[3, i] + rc[4, i])
        end
        out
    else
        throw(ArgumentError("mean!(out, f::BEDFile) requires out to be 1 by n or m by 1"))
    end
end

function missingrate!(out::AbstractMatrix{<:AbstractFloat}, f::BEDFile)
    k, l = size(out)
    m, n = size(f)
    if k == 1 && l == n
        cc = f.columncounts
        @inbounds for j in 1:n
            out[j] = cc[2, j] / m
        end
    elseif k == m && l == 1
        rc = _counts(f, 2)
        @inbounds for i in 1:m
            out[i] = rc[2, i] / n
        end
    else
        throw(ArgumentError("missingrate(out, f::BEDFile) requires out to be 1 by n or m by 1"))
    end
    out
end

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

function maffromcounts(v::SVector{4, Int})
    freq = (v[3] + 2v[4]) / 2(v[1] + v[3] + v[4])
    freq ≤ 0.5 ? freq : 1 - freq
end

function maf!(out::AbstractVector{T}, f::BEDFile) where T <: AbstractFloat
    cc = _counts(f, 1)
    map!(maffromcounts, out, f.staticcounts)
end

"""
    maf(f::BEDFile)

Return a vector of minor allele frequencies for the columns of `f`.

By definition the minor allele frequency is between 0 and 0.5
"""
maf(f::BEDFile) = maf!(Vector{Float64}(undef, size(f, 2)), f)

minorallelefromcounts(v::SVector{4, Int}) = v[1] > v[4]

function minorallele!(out::AbstractVector{Bool}, f::BEDFile)
    cc = _counts(f, 1)
    map!(minorallelefromcounts, out, f.staticcounts)
end

"""
    minorallele(f::BEDFile)

Return a `Vector{Bool}` indicating if the minor allele in each column is A2
"""
minorallele(f::BEDFile) = minorallele!(Vector{Bool}(undef, size(f, 2)), f)
