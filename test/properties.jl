
@testset "properties" begin

@testset "topology" begin
    cubane = smilestomol("C12C3C4C1C1C4C3C12")
    @test iterate(atom_sssrsizes(cubane)[1])[1] == 4
    @test iterate(bond_sssrsizes(cubane)[8])[1] == 4
    @test count(x->x==2, atom_sssrcount(cubane)) == 4
    @test count(x->x==2, bond_sssrcount(cubane)) == 8
    biphenyl = smilestomol("C1CCCCC1C1CCCCC1")
    @test isringatom(biphenyl)[6]
    @test isringatom(biphenyl)[7]
    @test !isringbond(biphenyl)[7]
end


@testset "elemental" begin
    alkyne = smilestomol("CC=CC#C")
    @test pielectron(alkyne) == [0, 1, 1, 2, 2]
    pyrrole = smilestomol("C1=CC=CN1")
    @test pielectron(pyrrole) == [1, 1, 1, 1, 0]
    azide = smilestomol("CC(=O)N=N=N")
    @test pielectron(azide) == [0, 1, 1, 1, 2, 1]
    amide = smilestomol("CCC(=O)N")
    @test ishacceptor(amide) == [0, 0, 0, 1, 1]
    @test ishdonor(amide) == [0, 0, 0, 0, 1]
    fluoro = smilestomol("CCN(CO)CF")
    @test ishacceptor(fluoro) == [0, 0, 1, 0, 1, 0, 1]
    @test ishdonor(fluoro) == [0, 0, 0, 0, 1, 0, 0]
end

@testset "rotatable" begin
    Phe = smilestomol("N[C@@H](CC1=CC=CC=C1)C(O)=O")
    @test rotatablecount(Phe) == 3
    KCl = smilestomol("[K+].[Cl-]")
    @test rotatablecount(KCl) == 0
    dipyridamole = smilestomol(
        "n3c(nc2c(nc(nc2N1CCCCC1)N(CCO)CCO)c3N4CCCCC4)N(CCO)CCO")
    @test rotatablecount(dipyridamole) == 12
    paclitaxel = smilestomol(
        "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
    )
    @test rotatablecount(paclitaxel) == 15
end

@testset "aromatic" begin
    phenol = smilestomol("C=1C=CC=CC=1O")  # not carbonyl O
    @test count(isaromatic(phenol)) == 6

    furan = smilestomol("o1cccc1")
    @test count(isaromatic(furan)) == 5

    quinone = smilestomol("C1(=O)C=CC(=O)C=C1")
    @test count(isaromatic(quinone)) == 0

    tropone = smilestomol("C1(=O)C=CC=CC=C1")
    @test count(isaromatic(tropone)) == 7

    azepine = smilestomol("N1C=CC=CC=C1")
    @test count(isaromatic(azepine)) == 0

    thiopheneoxide = smilestomol("C1=CC=CS1=O")
    @test count(isaromatic(thiopheneoxide)) == 0

    dihydroazaborine = smilestomol("C1=CC=CBN1")
    @test count(isaromatic(dihydroazaborine)) == 6

    fulvene = smilestomol("C1=CC=CC1=C")
    @test count(isaromatic(fulvene)) == 0

    cyclopropenyl = smilestomol("C1=C[C+]1")
    @test count(isaromatic(cyclopropenyl)) == 3

    hexatriene = smilestomol("C=CC=CC=C")
    @test count(isaromatic(hexatriene)) == 0

    borazine = smilestomol("B1NBNBN1")
    @test count(isaromatic(borazine)) == 6

    methylenedioxybenzene = smilestomol("C1=CC=CC2=C1OCO2")
    @test count(isaromatic(methylenedioxybenzene)) == 6

    naphthalene = smilestomol("C=1C=CC=C2C1C=CC=C2")
    @test count(isaromatic(naphthalene)) == 10

    coumarin = smilestomol("C1=CC(=O)OC2=C1C=CC=C2")
    @test count(isaromatic(coumarin)) == 10

    pyrene = smilestomol("C12=CC=C3C=CC=C4C=CC(C2=C34)=CC=C1")
    @test count(isaromatic(pyrene)) == 16

    caffeine = smilestomol("CN1C=NC2=C1C(=O)N(C)C(=O)N2C")
    @test count(isaromatic(caffeine)) == 9

    # Difficult aromaticities
    quinodimethane = smilestomol("C1=CC=CC(=C)C1=C")
    @test_broken count(isaromatic(azulene)) == 10

    dihydronaphthalene = smilestomol("C=1CCC=C2C1C=CC=C2")
    @test_broken count(isaromatic(dihydronaphthalene)) == 10

    azulene = smilestomol("C=1C=CC=2C1C=CC=CC2")
    @test_broken count(isaromatic(azulene)) == 10

    pyromellitimide = smilestomol("C1(=O)NC(=O)C=2C1=CC=3C(=O)NC(=O)C=3C=2")
    @test count(isaromatic(pyromellitimide)) == 6
end

end # properties
