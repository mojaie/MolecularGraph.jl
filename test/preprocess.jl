
@testset "preprocess" begin

@testset "allhydrogens" begin
    ethanol = smilestomol("[H]C([H])([H])C([H])([H])O")
    @test issetequal(allhydrogens(ethanol), [1, 3, 4, 6, 7])
    removed = removehydrogens(ethanol)
    @test nodecount(removed) == 3
    @test edgecount(removed) == 2
end

@testset "trivialhydrogens" begin
    ethanold3 = smilestomol("[2H]C([2H])([2H])C([H])([H])O")
    @test issetequal(trivialhydrogens(ethanold3), [6, 7])
    removed = removehydrogens(ethanold3, all=false)
    @test nodecount(removed) == 6
    @test edgecount(removed) == 5
    LAla = parse(SMILES, "[H]N([H])[C@][H](C)C(=O)O")
    @test issetequal(trivialhydrogens(LAla), [1, 3])
end

@testset "addhydrogens" begin
    neop = smilestomol("CC(C)(C)CO")
    addhydrogens!(neop)
    @test nodecount(neop) == 18
    @test edgecount(neop) == 17
end

@testset "largestcomponent" begin
    thiols = smilestomol("CCCCS.SCCCCC")
    @test issetequal(largestcomponentnodes(thiols), 6:11)
    lc = extractlargestcomponent(thiols)
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

@testset "kekulize" begin
    furan = smilestomol("o1cccc1")
    kekulize!(furan)
    @test edgeattr(furan, 2).order == 2
    @test edgeattr(furan, 4).order == 2

    pyrene = smilestomol("c1cc2cccc3ccc4cccc1c4c32")
    kekulize!(pyrene)
    @test sum(bondorder(pyrene) .== 1) == 11
    @test sum(bondorder(pyrene) .== 2) == 8
end

end # preprocess
