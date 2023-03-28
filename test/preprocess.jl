
@testset "preprocess" begin


@testset "kekulize" begin
    null = smilestomol("")
    @test isempty(bond_order(null))

    furan = smilestomol("o1cccc1")
    @test sum(bond_order(furan)) == 7

    pyrrole = smilestomol("[nH]1cccc1")
    @test sum(bond_order(pyrrole)) == 8
    @test_throws ErrorException smilestomol("n1cccc1")  # Wrong pyrrole

    pyrene = smilestomol("c1cc2cccc3ccc4cccc1c4c32")
    p = bond_order(pyrene)
    @test sum(p .== 1) == 11
    @test sum(p .== 2) == 8
    
    pyridineoxide = smilestomol("[n+]1([O-])ccccc1")
    @test sum(bond_order(pyridineoxide) .== 2) == 3

    c60 = smilestomol("c12c3c4c5c1c6c7c8c2c9c1c3c2c3c4c4c%10c5c5c6c6c7c7c%11c8c9c8c9c1c2c1c2c3c4c3c4c%10c5c5c6c6c7c7c%11c8c8c9c1c1c2c3c2c4c5c6c3c7c8c1c23")
    @test sum(bond_order(c60) .== 2) == 30
end

@testset "all_hydrogens" begin
    ethanol = smilestomol("[H]C([H])([H])C([H])([H])O")
    @test all_hydrogens(ethanol) == [1, 3, 4, 6, 7]
    remove_hydrogens!(ethanol)
    @test nv(ethanol) == 3
    @test ne(ethanol) == 2
end

@testset "removable_hydrogens" begin
    ethanold3 = smilestomol(
        MolGraph{Int,SMILESAtom,SMILESBond}, "[2H]C([2H])([2H])C([H])([H])O")
    @test removable_hydrogens(ethanold3) == [6, 7]
    remove_hydrogens!(ethanold3, all=false)
    @test nv(ethanold3) == 6
    @test ne(ethanold3) == 5
    LAla = smilestomol("[H]N([H])[C@][H](C)C(=O)O")
    @test removable_hydrogens(LAla) == [1, 3]
end

@testset "add_hydrogens" begin
    neop = smilestomol(MolGraph{Int,SMILESAtom,SMILESBond}, "CC(C)(C)CO")
    add_hydrogens!(neop)
    @test nv(neop) == 18
    @test ne(neop) == 17
end

@testset "largest_component" begin
    thiols = smilestomol(MolGraph{Int,SMILESAtom,SMILESBond}, "CCCCS.SCCCCC")
    @test largest_component_nodes(thiols) == collect(6:11)
    extract_largest_component!(thiols)
    @test nv(thiols) == 6
    @test ne(thiols) == 5
end

@testset "protonate_acids" begin
    AcOH = smilestomol("CC(=O)[O-]")
    @test protonate_acids(AcOH)[4] == 0

    thiol = smilestomol("CC[S-]")
    @test protonate_acids(thiol)[3] == 0

    noxide = smilestomol("[n+]1([O-])ccccc1")
    @test protonate_acids(noxide)[1] == 1
end

@testset "deprotonate_oniums" begin
    pyrrolidinium = smilestomol("C1[N+]CCC1")
    @test deprotonate_oniums(pyrrolidinium)[2] == 0

    oxonium = smilestomol("[H][O+]([H])[H]")
    @test deprotonate_oniums(oxonium)[2] == 0

    ammonium = smilestomol("C[N+](C)(C)C")
    @test deprotonate_oniums(ammonium)[2] == 1
end

@testset "polarize" begin
    acetone = smilestomol("C[C+]([O-])C")
    carr, oarr = depolarize(acetone)
    @test carr[[2, 3]] == [0, 0]
    @test oarr[2] == 2
    carr, oarr = polarize(acetone)
    @test carr[[2, 3]] == [1, -1]
    @test oarr[2] == 1

    phosphorylCl = smilestomol("[O-][P+](Cl)(Cl)Cl")
    carr, oarr = depolarize(phosphorylCl)
    @test carr[[1, 2]] == [0, 0]
    @test oarr[1] == 2
    carr, oarr = polarize(phosphorylCl)
    @test carr[[1, 2]] == [-1, 1]
    @test oarr[1] == 1

    noxide = smilestomol("[n]1(=O)ccccc1")  # incorrect n-oxide but valid
    carr, oarr = depolarize(noxide)
    @test carr[[1, 2]] == [0, 0]
    @test oarr[1] == 2
    carr, oarr = polarize(noxide)
    @test carr[[1, 2]] == [1, -1]
    @test oarr[1] == 1
    
    dmso = smilestomol("CS(=O)C")
    carr, oarr = polarize(dmso)
    @test carr[[2, 3]] == [1, -1]
    @test oarr[2] == 1
    carr, oarr = depolarize(dmso)
    @test carr[[2, 3]] == [0, 0]
    @test oarr[2] == 2
    carr, oarr = polarize(dmso, positive=Symbol[])
    @test carr[[2, 3]] == [0, 0]
    @test oarr[2] == 2
end

@testset "13dipole" begin
    me_azide = smilestomol("C[N-][N+]#N")
    carr, oarr = to_allene_like(me_azide)
    @test carr[[2, 4]] == [0, -1]
    @test oarr[[2, 3]] == [2, 2]
    carr, oarr = to_triple_bond(me_azide)
    @test carr[[2, 4]] == [-1, 0]
    @test oarr[[2, 3]] == [1, 3]
end

end # preprocess
