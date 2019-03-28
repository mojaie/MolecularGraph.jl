
@testset "preprocess" begin

@testset "all_hydrogens" begin
    ethanol = smilestomol("[H]C([H])([H])C([H])([H])O")
    @test issetequal(all_hydrogens(ethanol), [1, 3, 4, 6, 7])
    removed = make_hydrogens_implicit(ethanol)
    @test nodecount(removed) == 3
    @test edgecount(removed) == 2
end

@testset "trivial_hydrogens" begin
    ethanol = smilestomol("[H]C([H])([H])C([H])([H])O")
    @test issetequal(trivialhydrogens(ethanol), [1, 3, 4, 6, 7])
    removed = make_hydrogens_implicit(ethanol, all=false)
    @test nodecount(removed) == 3
    @test edgecount(removed) == 2
end

@testset "hydrogens_explicit" begin
    neop = smilestomol("CC(C)(C)CO")
    added = make_hydrogens_explicit(neop)
    @test nodecount(added) == 18
    @test edgecount(added) == 17
end

@testset "largest_component" begin
    thiols = smilestomol("CCCCS.SCCCCC")
    @test issetequal(largest_component_nodes(thiols), 6:11)
    lc = largestcomponent(thiols)
    @test nodecount(lc) == 6
    @test edgecount(lc) == 5
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

@testset "depolarize" begin
    acetone = smilestomol("C[C+]([O-])C")
    depolarize!(acetone)
    @test getnode(acetone, 2).charge == 0
    @test getnode(acetone, 3).charge == 0
    @test getedge(acetone, 2).order == 2
end

@testset "triplebond_anion" begin
    me_azide = smilestomol("C[N-][N+]#N")
    triplebond_anion!(me_azide)
    @test getnode(me_azide, 2).charge == 0
    @test getnode(me_azide, 4).charge == -1
    @test getedge(me_azide, 2).order == 2
    @test getedge(me_azide, 3).order == 2
end

end # preprocess
