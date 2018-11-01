using BEDFiles, SparseArrays, Test

const EUR = BEDFile(BEDFiles.datadir("EUR_subset.bed"))
const mouse = BEDFile(BEDFiles.datadir("mouse.bed"))

@testset "size" begin
    @test size(EUR) == (379, 54051)
    @test size(EUR, 1) == 379
    @test size(EUR, 2) == 54051
    @test size(mouse) == (1940, 10150)
    @test size(mouse, 3) == 1
    @test_throws ErrorException size(mouse, 0)
    @test length(mouse) == 1940 * 10150
    @test eltype(mouse) == UInt8
end

@testset "missingpos" begin
    mp = missingpos(mouse)
    cc = counts(mouse, dims=1)
    @test isa(mp, SparseMatrixCSC)
    @test sum(mp, dims = 1) == view(cc, 2:2, :)
    @test sum(mp, dims = 2) == view(counts(mouse, dims=2), 2:2, :)'
end

include("summaries.jl")
