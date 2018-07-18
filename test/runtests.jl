using BEDFiles, Test

const EUR = BEDFile(BEDFiles.datadir("EUR_subset.bed"))

@testset "size" begin
    @test size(EUR) == (379, 54051)
    @test size(EUR, 1) == 379
    @test size(EUR, 2) == 54051
end

@testset "counts" begin
    cc = counts(EUR, dims=1)
    @test size(cc) == (4, size(EUR, 2))
    @test view(cc, :, 1) == [2, 0, 70, 307]
    @test all(sum(cc, dims = 1) .== size(EUR, 1))
    @test all(iszero, view(cc, 2, :))
    rc = counts(EUR, dims=2)
    @test size(rc) == (4, size(EUR, 1))
    @test view(rc, :, 1) == [2997, 0, 13143, 37911]
    @test all(sum(rc, dims = 1) .== size(EUR, 2))
    @test all(iszero, view(rc, 2, :))
end

@testset "means" begin
    cmns = mean(EUR, dims = 1)
    @test size(cmns) == (1, size(EUR, 2))
    @test cmns[1] ≈ 1.804749340369393
    @test mean(cmns) ≈ mean(EUR)
    @test minimum(cmns) ≈ 1.0
    @test maximum(cmns) ≈ 1.9788918205804749
    rmns = mean(EUR, dims = 2)
    @test size(rmns) == (size(EUR, 1), 1)
    @test rmns[1] ≈ 1.6459454959205195
    @test mean(rmns) ≈ mean(EUR)
    rmnextrma = extrema(rmns)
    @test rmnextrma[1] ≈ 1.637749532848606
    @test rmnextrma[2] ≈ 1.6627259440158368
end
