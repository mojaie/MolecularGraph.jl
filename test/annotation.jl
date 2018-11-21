
@testset "annotation" begin

@testset "valence" begin
    alkyne = smilestomol("CC=CC#C")
    @test alkyne.v[:Pi] == [0, 1, 1, 2, 2]
    pyrrole = smilestomol("C1=CC=CN1")
    @test pyrrole.v[:Pi] == [1, 1, 1, 1, 0]
    azide = smilestomol("CC(=O)N=N=N")
    @test azide.v[:Pi] == [0, 1, 1, 1, 2, 1]
    amide = smilestomol("CCC(=O)N")
    @test amide.v[:H_Acceptor] == [0, 0, 0, 1, 1]
    @test amide.v[:H_Donor] == [0, 0, 0, 0, 1]
    fluoro = smilestomol("CCN(CO)CF")
    @test fluoro.v[:H_Acceptor] == [0, 0, 1, 0, 1, 0, 1]
    @test fluoro.v[:H_Donor] == [0, 0, 0, 0, 1, 0, 0]
end

@testset "group_annot" begin
    assign_annot!(mol) = (atom_annot!(mol); group_annot!(mol))
    nullmol = smilestomol("")
    assign_annot!(nullmol)
    @test count(nullmol.v[:Carbonyl]) == 0
    alcohol = smilestomol("CC(O)C(CO)(C)O")
    assign_annot!(alcohol)
    @test alcohol.v[:Oxygen] == [
        nothing, nothing, 1, nothing, nothing, 1, nothing, 1]
    @test alcohol.v[:Alcohol] == [
        nothing, 2, nothing, 3, 1, nothing, nothing, nothing]
    amine = smilestomol("CNC(C#N)(C=N)N")
    assign_annot!(amine)
    @test amine.v[:Nitrogen] == [
        nothing, 2, nothing, nothing, 8, nothing, 5, 1]
end

@testset "rotatable" begin
    Phe = smilestomol("N[C@@H](CC1=CC=CC=C1)C(O)=O")
    rotatable!(Phe)
    @test count(Phe.v[:Rotatable]) == 3
    KCl = smilestomol("[K+].[Cl-]")
    rotatable!(KCl)
    @test count(KCl.v[:Rotatable]) == 0
    dipyridamole = smilestomol(
        "n3c(nc2c(nc(nc2N1CCCCC1)N(CCO)CCO)c3N4CCCCC4)N(CCO)CCO")
    rotatable!(dipyridamole)
    @test count(dipyridamole.v[:Rotatable]) == 12
    paclitaxel = smilestomol(
        "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
    )
    rotatable!(paclitaxel)
    @test count(paclitaxel.v[:Rotatable]) == 15
end

@testset "aromatic" begin
    assign_annot!(mol) = (atom_annot!(mol); group_annot!(mol); aromatic!(mol))
    furan = smilestomol("O1C=CC=C1")
    assign_annot!(furan)
    @test count(furan.v[:Aromatic]) == 5
    quinone = smilestomol("C1(=O)C=CC(=O)C=C1")
    assign_annot!(quinone)
    @test count(quinone.v[:Aromatic]) == 0
    tropone = smilestomol("C1(=O)C=CC=CC=C1")
    assign_annot!(tropone)
    @test count(tropone.v[:Aromatic]) == 7
    azepine = smilestomol("N1C=CC=CC=C1")
    assign_annot!(azepine)
    @test count(azepine.v[:Aromatic]) == 0
    thiopheneoxide = smilestomol("C1=CC=CS1=O")
    assign_annot!(thiopheneoxide)
    @test count(thiopheneoxide.v[:Aromatic]) == 0
    dihydroazaborine = smilestomol("C1=CC=CBN1")
    assign_annot!(dihydroazaborine)
    @test count(dihydroazaborine.v[:Aromatic]) == 6
    fulvene = smilestomol("C1=CC=CC1=C")
    assign_annot!(fulvene)
    @test count(fulvene.v[:Aromatic]) == 0
    cyclopropenyl = smilestomol("C1=C[C+]1")
    assign_annot!(cyclopropenyl)
    @test count(cyclopropenyl.v[:Aromatic]) == 3
    hexatriene = smilestomol("C=CC=CC=C")
    assign_annot!(hexatriene)
    @test count(hexatriene.v[:Aromatic]) == 0
    borazine = smilestomol("B1NBNBN1")
    assign_annot!(borazine)
    @test count(borazine.v[:Aromatic]) == 6
    methylenedioxybenzene = smilestomol("C1=CC=CC2=C1OCO2")
    assign_annot!(methylenedioxybenzene)
    @test count(methylenedioxybenzene.v[:Aromatic]) == 6
    naphthalene = smilestomol("C=1C=CC=C2C1C=CC=C2")
    assign_annot!(naphthalene)
    @test count(naphthalene.v[:Aromatic]) == 10
    coumarin = smilestomol("C1=CC(=O)OC2=C1C=CC=C2")
    assign_annot!(coumarin)
    @test count(coumarin.v[:Aromatic]) == 10
    pyrene = smilestomol("C12=CC=C3C=CC=C4C=CC(C2=C34)=CC=C1")
    assign_annot!(pyrene)
    @test count(pyrene.v[:Aromatic]) == 16
    caffeine = smilestomol("CN1C=NC2=C1C(=O)N(C)C(=O)N2C")
    assign_annot!(caffeine)
    @test count(caffeine.v[:Aromatic]) == 9

    # Difficult aromaticities
    quinodimethane = smilestomol("C1=CC=CC(=C)C1=C")
    assign_annot!(quinodimethane)
    @test_broken count(azulene.v[:Aromatic]) == 10
    dihydronaphthalene = smilestomol("C=1CCC=C2C1C=CC=C2")
    assign_annot!(dihydronaphthalene)
    @test_broken count(dihydronaphthalene.v[:Aromatic]) == 10
    azulene = smilestomol("C=1C=CC=2C1C=CC=CC2")
    assign_annot!(azulene)
    @test_broken count(azulene.v[:Aromatic]) == 10
    pyromellitimide = smilestomol("C1(=O)NC(=O)C=2C1=CC=3C(=O)NC(=O)C=3C=2")
    assign_annot!(pyromellitimide)
    @test count(pyromellitimide.v[:Aromatic]) == 6
end

end # basedescriptor
