
@testset "preprocess" begin

@testset "kekulize" begin
    furan = smilestomol("o1cccc1")
    furan = kekulize(furan)
    @test edgeattr(furan, 2).order == 2
    @test edgeattr(furan, 4).order == 2

    pyrene = smilestomol("c1cc2cccc3ccc4cccc1c4c32")
    pyrene = kekulize(pyrene)
    @test sum(bondorder(pyrene) .== 1) == 11
    @test sum(bondorder(pyrene) .== 2) == 8
    
    pyridineoxide = smilestomol("[n+]1([O-])ccccc1")
    pyridineoxide = kekulize(pyridineoxide)
    @test sum(bondorder(pyridineoxide) .== 2) == 3
end

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
    neop = addhydrogens(neop)
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

@testset "protonateacids" begin
    AcOH = smilestomol("CC(=O)[O-]")
    AcOH = protonateacids(AcOH)
    @test nodeattr(AcOH, 4).charge == 0

    thiol = smilestomol("CC[S-]")
    thiol = protonateacids(thiol)
    @test nodeattr(thiol, 3).charge == 0

    noxide = smilestomol("[n+]1([O-])ccccc1")
    noxide = protonateacids(noxide)
    @test nodeattr(noxide, 1).charge == 1
end

@testset "deprotonateoniums" begin
    pyrrolidinium = smilestomol("C1[N+]CCC1")
    pyrrolidinium = deprotonateoniums(pyrrolidinium)
    @test nodeattr(pyrrolidinium, 2).charge == 0

    oxonium = smilestomol("[H][O+]([H])[H]")
    oxonium = deprotonateoniums(oxonium)
    @test nodeattr(oxonium, 2).charge == 0

    ammonium = smilestomol("C[N+](C)(C)C")
    ammonium = deprotonateoniums(ammonium)
    @test nodeattr(ammonium, 2).charge == 1
end

@testset "polarize" begin
    acetone = smilestomol("C[C+]([O-])C")
    acetone = depolarize(acetone)
    @test nodeattr(acetone, 2).charge == 0
    @test nodeattr(acetone, 3).charge == 0
    @test edgeattr(acetone, 2).order == 2
    acetone = polarize(acetone)
    @test nodeattr(acetone, 2).charge == 0
    @test nodeattr(acetone, 3).charge == 0
    @test edgeattr(acetone, 2).order == 2

    phosphorylCl = smilestomol("[O-][P+](Cl)(Cl)Cl")
    phosphorylCl = depolarize(phosphorylCl)
    @test nodeattr(phosphorylCl, 1).charge == 0
    @test nodeattr(phosphorylCl, 1).charge == 0
    @test edgeattr(phosphorylCl, 1).order == 2
    phosphorylCl = polarize(phosphorylCl)
    @test nodeattr(phosphorylCl, 1).charge == 0
    @test nodeattr(phosphorylCl, 1).charge == 0
    @test edgeattr(phosphorylCl, 1).order == 2

    noxide = smilestomol("[n]1(=O)ccccc1")
    noxide = polarize(noxide)
    @test nodeattr(noxide, 1).charge == 1
    @test nodeattr(noxide, 2).charge == -1
    @test edgeattr(noxide, 1).order == 1
    noxide = depolarize(noxide)
    @test nodeattr(noxide, 1).charge == 1
    @test nodeattr(noxide, 2).charge == -1
    @test edgeattr(noxide, 1).order == 1
    
    dmso = smilestomol("CS(=O)C")
    dmso = polarize(dmso)
    @test nodeattr(dmso, 2).charge == 1
    @test nodeattr(dmso, 3).charge == -1
    @test edgeattr(dmso, 2).order == 1
    dmso = depolarize(dmso)
    @test nodeattr(dmso, 2).charge == 1
    @test nodeattr(dmso, 3).charge == -1
    @test edgeattr(dmso, 2).order == 1
    dmso = depolarize(dmso, positive=[:S])
    @test nodeattr(dmso, 2).charge == 0
    @test nodeattr(dmso, 3).charge == 0
    @test edgeattr(dmso, 2).order == 2
end

@testset "13dipole" begin
    me_azide = smilestomol("C[N-][N+]#N")
    me_azide = toallenelike(me_azide)
    @test nodeattr(me_azide, 2).charge == 0
    @test nodeattr(me_azide, 4).charge == -1
    @test edgeattr(me_azide, 2).order == 2
    @test edgeattr(me_azide, 3).order == 2
    me_azide = totriplebond(me_azide)
    @test nodeattr(me_azide, 2).charge == -1
    @test nodeattr(me_azide, 4).charge == 0
    @test edgeattr(me_azide, 2).order == 1
    @test edgeattr(me_azide, 3).order == 3
end

end # preprocess
