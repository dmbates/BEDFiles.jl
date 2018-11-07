"""
    uvals(v::SVector{4,Int}; impute::Bool=true, center::Bool=true, scale::Bool=true)

Return an `SVector{4,Float64}` of the unique values in the additive model incidence
"""
function uvals!(uv::Vector{<:AbstractFloat}, v::AbstractVector{<:Integer})
    μ = meanfromcounts(v)
    σ = sqrt(μ * (1 - μ / 2))
    uv[1] = -μ / σ
    uv[2] = 0
    uv[3] = (1 - μ) / σ
    uv[4] = (2 - μ) / σ
    return uv
end

function newgrm(f::BEDFile, tempm::AbstractMatrix{<:AbstractFloat};
        colinds::Union{Nothing,AbstractVector{Bool}}=nothing,
        mafthreshold::Real=0.01)
    m, n = size(f)
    k, l = size(tempm)
    m == k || throw(DimensionMismatch("size(f, 1) ≠ size(tempm, 1)"))
    T = eltype(tempm)
    uvec = Vector{T}(undef, 4)
    result = zeros(T, (m, m))
    colinds = something(colinds, maffromcounts.(f.sccounts) .≥ mafthreshold)
    sccounts = f.sccounts
    ind = 1
    while ind ≤ n
        j = 1
        while j ≤ l && ind ≤ n
            if colinds[ind]
                uvals!(uvec, sccounts[ind])
                @inbounds for i in 1:m
                    tempm[i, j] = uvec[f[i,ind] + 1]
                end
                j += 1
            end
            ind += 1
        end
        BLAS.syrk!('L', 'N', one(T), j < l ? view(tempm, :, 1:j) : tempm, one(T), result)
    end
    result ./= 2n
    Symmetric(result, :L)
end
