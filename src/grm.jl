"""
    uvals(v::SVector{4,Int}; impute::Bool=true, center::Bool=true, scale::Bool=true)

Return an `SVector{4,Float64}` of the unique values in the additive model incidence
"""
function uvals!(uv::Vector{<:AbstractFloat}, v::SVector{4,Int}; impute::Bool=true, center::Bool=true, scale::Bool=true)
    μ = meanfromcounts(v)
    σ = sqrt(μ * (1 - μ / 2))
    uv[1] = -μ / σ
    uv[2] = 0
    uv[3] = (1 - μ) / σ
    uv[4] = (2 - μ) / σ
    uv
end

function newgrm(f::BEDFile)
    m, n = size(f)
    uvec = Vector{Float64}(undef, 4)
    tempm = Matrix{Float64}(undef, size(f))
    @inbounds for (j, v) in enumerate(f.sccounts)
        uvals!(uvec, v)
        for i in 1:m
            tempm[i, j] = uvec[f[i,j] + 1]
        end
    end
    (tempm * tempm') ./ (2n)
end