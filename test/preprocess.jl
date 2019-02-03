
@testset "preprocess" begin

@testset "allhydrogens" begin
    ethanol = parse(SMILES, "[H]C([H])([H])C([H])([H])O")
    @test issetequal(allhydrogens(ethanol), [1, 3, 4, 6, 7])
end

@testset "trivialhydrogens" begin
    ethanol = parse(SMILES, "[H]C([H])([H])C([H])([H])O")
    @test issetequal(trivialhydrogens(ethanol), [1, 3, 4, 6, 7])
end

@testset "largestcomponent" begin
    thiols = parse(SMILES, "CCCCS.SCCCCC")
    @test issetequal(largestcomponent(thiols), 6:11)
end

@testset "neutralize_acids" begin
    AcOH = smilestomol("CC(=O)[O-]")
    neutralize_acids!(AcOH)
    @test getnode(AcOH, 4).charge == 0
end

@testset "neutralize_oniums" begin
    pyrrolidinium = smilestomol("C1[N+]CCC1")
    neutralize_oniums!(pyrrolidinium)
    @test getnode(pyrrolidinium, 2).charge == 0
end

end # preprocess
