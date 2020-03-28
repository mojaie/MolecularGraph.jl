
@testset "preprocess" begin

@testset "all_hydrogens" begin
    ethanol = smilestomol("[H]C([H])([H])C([H])([H])O")
    @test issetequal(allhydrogens(ethanol), [1, 3, 4, 6, 7])
    removed = removehydrogens(ethanol)
    @test nodecount(removed) == 3
    @test edgecount(removed) == 2
end

@testset "trivial_hydrogens" begin
    # TODO:
    ethanol = smilestomol("[H]C([H])([H])C([H])([H])O")
    @test issetequal(trivialhydrogens(ethanol), [1, 3, 4, 6, 7])
    removed = removehydrogens(ethanol, all=false)
    @test nodecount(removed) == 3
    @test edgecount(removed) == 2
end

@testset "addhydrogens" begin
    neop = smilestomol("CC(C)(C)CO")
    addhydrogens!(neop)
    @test nodecount(neop) == 18
    @test edgecount(neop) == 17
end

@testset "largest_component" begin
    thiols = smilestomol("CCCCS.SCCCCC")
    @test issetequal(largestcomponentnodes(thiols), 6:11)
    lc = largestcomponentgraph(thiols)
    @test nodecount(lc) == 6
    @test edgecount(lc) == 5
end

@testset "neutralize_acids" begin
    AcOH = smilestomol("CC(=O)[O-]")
    neutralizeacids!(AcOH)
    @test nodeattr(AcOH, 4).charge == 0
end

@testset "neutralize_oniums" begin
    pyrrolidinium = smilestomol("C1[N+]CCC1")
    neutralizeoniums!(pyrrolidinium)
    @test nodeattr(pyrrolidinium, 2).charge == 0
end

@testset "depolarize" begin
    acetone = smilestomol("C[C+]([O-])C")
    depolarize!(acetone)
    @test nodeattr(acetone, 2).charge == 0
    @test nodeattr(acetone, 3).charge == 0
    @test edgeattr(acetone, 2).order == 2
end

@testset "toallenelike" begin
    me_azide = smilestomol("C[N-][N+]#N")
    toallenelike!(me_azide)
    @test nodeattr(me_azide, 2).charge == 0
    @test nodeattr(me_azide, 4).charge == -1
    @test edgeattr(me_azide, 2).order == 2
    @test edgeattr(me_azide, 3).order == 2
end

end # preprocess
