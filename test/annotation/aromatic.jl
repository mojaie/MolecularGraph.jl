
@testset "annotation.aromatic" begin

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

end # annotation.aromatic
