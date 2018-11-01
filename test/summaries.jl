using BEDFiles, Test

if !@isdefined(EUR) || !isa(EUR, BEDFile)
    const EUR = BEDFile(BEDFiles.datadir("EUR_subset.bed"))
end;
if !@isdefined(mouse) || !isa(mouse, BEDFile)
    const mouse = BEDFile(BEDFiles.datadir("mouse.bed"))
end;

@testset "counts" begin
    cc = counts(EUR, dims=1)
    @test counts(EUR) == vec(sum(cc, dims=2))
    @test size(cc) == (4, size(EUR, 2))
    @test view(cc, :, 1) == [2, 0, 70, 307]
    @test all(sum(cc, dims = 1) .== size(EUR, 1))
    @test all(iszero, view(cc, 2, :))
    rc = counts(EUR, dims=2)
    @test size(rc) == (4, size(EUR, 1))
    @test view(rc, :, 1) == [2997, 0, 13143, 37911]
    @test all(sum(rc, dims = 1) .== size(EUR, 2))
    @test all(iszero, view(rc, 2, :))
    @test_throws ArgumentError counts(EUR, dims=3)
    @test sum(counts(EUR)) == length(EUR)
end

@testset "means" begin
    cmns = mean(EUR, dims = 1)
    @test size(cmns) == (1, size(EUR, 2))
    @test cmns ≈ mean!(similar(cmns), EUR)
    @test cmns[1] ≈ 1.804749340369393
    @test mean(EUR) ≈ 1.6497440680596342
    @test minimum(cmns) ≈ 1.0
    @test maximum(cmns) ≈ 1.9788918205804749
    rmns = mean(EUR, dims = 2)
    @test size(rmns) == (size(EUR, 1), 1)
    @test rmns[1] ≈ 1.6459454959205195
    @test rmns ≈ mean!(Matrix{Float64}(undef, size(rmns)), EUR)
    @test mean(rmns) ≈ mean(EUR)
    rmnextrma = extrema(rmns)
    @test rmnextrma[1] ≈ 1.637749532848606
    @test rmnextrma[2] ≈ 1.6627259440158368
    @test_throws ArgumentError mean(EUR, dims=3)
    @test mean(EUR) ≈ mean!(Matrix{Float64}(undef, (1,1)), EUR)[1]
end

@testset "missingrate" begin
    @test iszero(missingrate(EUR))
    EURmisscols = missingrate(EUR, dims=1)
    @test size(EURmisscols) == (1, size(EUR, 2))
    @test all(iszero, EURmisscols)
    @test all(iszero, missingrate!(similar(EURmisscols), EUR))
    EURmissrows = missingrate(EUR, dims=2)
    @test size(EURmissrows) == (size(EUR, 1), 1)
    @test all(iszero, EURmissrows)
    @test EURmissrows == missingrate!(similar(EURmissrows), EUR)
    @test missingrate(mouse) ≈ 0.0017227159616068255
    @test missingrate!(Matrix{Float64}(undef, (1,1)), mouse)[1] ≈ 0.0017227159616068255
    mousemcols = missingrate(mouse, dims=1)
    @test size(mousemcols) == (1, size(mouse, 2))
    @test mousemcols ≈ missingrate!(similar(mousemcols), mouse)
    @test iszero(minimum(mousemcols))
    @test maximum(mousemcols) ≈ 0.09845360824742268
end

@testset "var/std" begin
    cvars = var(EUR, dims = 1)
    @test size(cvars) == (1, size(EUR, 2))
    @test cvars[1] ≈ 0.16812553224162724
    @test iszero(minimum(cvars))
    @test maximum(cvars) ≈ 0.8705727966941688
    @test std(EUR, dims = 1) ≈ sqrt.(cvars)
    rvars = var(EUR, dims = 2)
    @test size(rvars) == (size(EUR, 1), 1)
    @test rvars[1] ≈ 0.33960146078503206
    @test minimum(rvars) ≈ 0.31972707244268267
    @test maximum(rvars) ≈ 0.365901927927013
    @test std(EUR, dims = 2) ≈ sqrt.(rvars)
end

@testset "minorallele" begin
    mall = minorallele(EUR)
    @test !any(mall)
    @test !any(minorallele!(similar(mall), EUR))
    mall = minorallele(mouse)
    @test !any(mall)
    @test !any(minorallele!(similar(mall), mouse))
end

@testset "maf" begin
    mafv = maf(EUR)
    @test size(mafv) == (size(EUR, 2),)
    @test mafv ≈ maf!(similar(mafv), EUR)
    @test minimum(mafv) ≈ 0.01055408970976257
    @test maximum(mafv) == 0.5
    @test mafv[1] ≈ 0.09762532981530347
    @test mafv[end] ≈ 0.02638522427440637
end