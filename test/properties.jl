
@testset "properties" begin

"""
@testset "topology" begin
    RESOURCE_DIR = joinpath(dirname(@__FILE__), "..", "_resources", "DrugBank")
    phe = loadsdfmol(open(joinpath(RESOURCE_DIR, "Phe.mol")))
    @test phe.rings[1].arr == UInt16[6, 7, 10, 12, 11, 8]
    @test phe.scaffolds == [[1]]
    @test length(phe.isolated) == 0
    prem = loadsdfmol(open(joinpath(RESOURCE_DIR, "Premarin.mol")))
    @test prem.scaffolds == [[1, 2, 3, 4]]
    @test prem.isolated == [[28]]
end
"""

@testset "elemental" begin
    alkyne = smilestomol("CC=CC#C")
    @test alkyne[:pielectron] == [0, 1, 1, 2, 2]
    pyrrole = smilestomol("C1=CC=CN1")
    @test pyrrole[:pielectron] == [1, 1, 1, 1, 0]
    azide = smilestomol("CC(=O)N=N=N")
    @test azide[:pielectron] == [0, 1, 1, 1, 2, 1]
    amide = smilestomol("CCC(=O)N")
    @test amide[:ishacceptor] == [0, 0, 0, 1, 1]
    @test amide[:ishdonor] == [0, 0, 0, 0, 1]
    fluoro = smilestomol("CCN(CO)CF")
    @test fluoro[:ishacceptor] == [0, 0, 1, 0, 1, 0, 1]
    @test fluoro[:ishdonor] == [0, 0, 0, 0, 1, 0, 0]
end

@testset "rotatable" begin
    Phe = smilestomol("N[C@@H](CC1=CC=CC=C1)C(O)=O")
    @test count(Phe[:isrotatable]) == 3
    KCl = smilestomol("[K+].[Cl-]")
    @test count(KCl[:isrotatable]) == 0
    dipyridamole = smilestomol(
        "n3c(nc2c(nc(nc2N1CCCCC1)N(CCO)CCO)c3N4CCCCC4)N(CCO)CCO")
    @test count(dipyridamole[:isrotatable]) == 12
    paclitaxel = smilestomol(
        "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
    )
    @test count(paclitaxel[:isrotatable]) == 15
end

@testset "aromatic" begin
    phenol = smilestomol("C=1C=CC=CC=1O")  # not carbonyl O
    @test count(phenol[:isaromatic]) == 6

    furan = smilestomol("o1cccc1")
    @test count(furan[:isaromatic]) == 5

    quinone = smilestomol("C1(=O)C=CC(=O)C=C1")
    @test count(quinone[:isaromatic]) == 0

    tropone = smilestomol("C1(=O)C=CC=CC=C1")
    @test count(tropone[:isaromatic]) == 7

    azepine = smilestomol("N1C=CC=CC=C1")
    @test count(azepine[:isaromatic]) == 0

    thiopheneoxide = smilestomol("C1=CC=CS1=O")
    @test count(thiopheneoxide[:isaromatic]) == 0

    dihydroazaborine = smilestomol("C1=CC=CBN1")
    @test count(dihydroazaborine[:isaromatic]) == 6

    fulvene = smilestomol("C1=CC=CC1=C")
    @test count(fulvene[:isaromatic]) == 0

    cyclopropenyl = smilestomol("C1=C[C+]1")
    @test count(cyclopropenyl[:isaromatic]) == 3

    hexatriene = smilestomol("C=CC=CC=C")
    @test count(hexatriene[:isaromatic]) == 0

    borazine = smilestomol("B1NBNBN1")
    @test count(borazine[:isaromatic]) == 6

    methylenedioxybenzene = smilestomol("C1=CC=CC2=C1OCO2")
    @test count(methylenedioxybenzene[:isaromatic]) == 6

    naphthalene = smilestomol("C=1C=CC=C2C1C=CC=C2")
    @test count(naphthalene[:isaromatic]) == 10

    coumarin = smilestomol("C1=CC(=O)OC2=C1C=CC=C2")
    @test count(coumarin[:isaromatic]) == 10

    pyrene = smilestomol("C12=CC=C3C=CC=C4C=CC(C2=C34)=CC=C1")
    @test count(pyrene[:isaromatic]) == 16

    caffeine = smilestomol("CN1C=NC2=C1C(=O)N(C)C(=O)N2C")
    @test count(caffeine[:isaromatic]) == 9

    # Difficult aromaticities
    quinodimethane = smilestomol("C1=CC=CC(=C)C1=C")
    @test_broken count(azulene[:isaromatic]) == 10

    dihydronaphthalene = smilestomol("C=1CCC=C2C1C=CC=C2")
    @test_broken count(dihydronaphthalene[:isaromatic]) == 10

    azulene = smilestomol("C=1C=CC=2C1C=CC=CC2")
    @test_broken count(azulene[:isaromatic]) == 10

    pyromellitimide = smilestomol("C1(=O)NC(=O)C=2C1=CC=3C(=O)NC(=O)C=3C=2")
    @test count(pyromellitimide[:isaromatic]) == 6
end

end # properties
