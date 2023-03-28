
@testset "properties" begin

@testset "topology" begin
    cubane = smilestomol("C12C3C4C1C1C4C3C12")
    @test all(is_in_ring(cubane))
    @test smallest_ring(cubane)[4] == 4
    biphenyl = smilestomol("C1CCCCC1C1CCCCC1")
    @test ring_count(biphenyl)[12] == 1
    @test is_in_ring(biphenyl)[6]
    @test is_in_ring(biphenyl)[7]
    @test !is_edge_in_ring(biphenyl)[edge_rank(biphenyl, 6, 7)]
end

@testset "elements" begin
    sodiumsulfate = smilestomol("O=S(=O)([O-])[O-].[Na+].[Na+]")
    @test atom_symbol(sodiumsulfate) == [:O, :S, :O, :O, :O, :Na, :Na]
    @test bond_order(sodiumsulfate) == [2, 2, 1, 1]
    @test charge(sodiumsulfate) == [0, 0, 0, -1, -1, 1, 1]
    @test degree(sodiumsulfate) == [1, 4, 1, 1, 1, 0, 0]

    isocyanurate = smilestomol("ClN1C(=O)N(Cl)C(=O)N(Cl)C1=O")
    @test atom_symbol(isocyanurate) == [
        :Cl, :N, :C, :O, :N, :Cl, :C, :O, :N, :Cl, :C, :O]
    @test bond_order(isocyanurate) == [
        1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2]
    @test all(charge(isocyanurate) .== 0)
    @test degree(isocyanurate) == [1, 3, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1]
end

@testset "formula" begin
    hexose = smilestomol("OCC(O)C(O)C(O)C(O)C=O")
    @test molecular_formula(hexose) == "C6H12O6"
    @test empirical_formula(hexose) == "CH2O"

    sulfalicacid = smilestomol("O=S(=O)(O)O")
    @test molecular_formula(sulfalicacid) == "H2O4S"
    @test empirical_formula(sulfalicacid) == "H2O4S"

    order = smilestomol("[Na].[Ba].B.Br.B.[Ar].Br")
    @test molecular_formula(order) == "ArB2BaBr2H8Na"
end

@testset "valence" begin
    atoms = smilestomol("B.C.N.O.F.[Si].P.S.Cl.[As].[Se].Br.I")
    @test valence(atoms) == [3, 4, 3, 2, 1, 4, 3, 2, 1, 3, 2, 1, 1]
    @test lone_pair(atoms) == [-1, 0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3, 3]
    @test implicit_hydrogens(atoms) == [3, 4, 3, 2, 1, 4, 3, 2, 1, 3, 2, 1, 1]
    
    pyridineoxide = smilestomol("[N+]1([O-])=CC=CC=C1")
    @test valence(pyridineoxide) == [4, 1, 4, 4, 4, 4, 4]
    @test lone_pair(pyridineoxide) == [0, 3, 0, 0, 0, 0, 0]
    @test implicit_hydrogens(pyridineoxide) == [0, 0, 1, 1, 1, 1, 1]
    
    trifluoroborate = smilestomol("C[B-](F)(F)F")
    @test valence(trifluoroborate) == [4, 4, 1, 1, 1]
    @test lone_pair(trifluoroborate) == [0, 0, 3, 3, 3]
    @test implicit_hydrogens(trifluoroborate) == [3, 0, 0, 0, 0]

    etmacl = smilestomol("[H]C([H])([H])C[Mg][Cl]")
    @test valence(etmacl) == [1, 4, 1, 1, 4, nothing, 1]
    @test lone_pair(etmacl) == [0, 0, 0, 0, 0, nothing, 3]
    @test implicit_hydrogens(etmacl) == [0, 0, 0, 0, 2, 0, 0]
    @test explicit_hydrogens(etmacl) == [0, 3, 0, 0, 0, 0, 0]
    @test total_hydrogens(etmacl) == [0, 3, 0, 0, 2, 0, 0]
    @test heavy_atoms(etmacl) == [1, 1, 1, 1, 2, 2, 1]
    @test connectivity(etmacl) == [1, 4, 1, 1, 4, 2, 1]
end

@testset "hybridization" begin
    alkyne = smilestomol("CC=CC#C")
    @test pi_electron(alkyne) == [0, 1, 1, 2, 2]
    @test hybridization(alkyne) == [:sp3, :sp2, :sp2, :sp, :sp]
    
    anilinium = smilestomol("C1=CC=CC=C1[N+]")
    @test pi_electron(anilinium) == [1, 1, 1, 1, 1, 1, 0]
    @test hybridization(anilinium) == [:sp2, :sp2, :sp2, :sp2, :sp2, :sp2, :sp3]

    azide = smilestomol("CC(=O)N=[N+]=[N-]")
    @test pi_electron(azide) == [0, 1, 1, 1, 2, 1]
    @test hybridization(azide) == [:sp3, :sp2, :sp2, :sp2, :sp, :sp2]

    pyrrole = smilestomol("C1=CC=CN1")
    @test pi_electron(pyrrole) == [1, 1, 1, 1, 2]
    @test hybridization(pyrrole) == [:sp2, :sp2, :sp2, :sp2, :sp2]

    nacl = smilestomol("[Na+][Cl-]")
    @test pi_electron(nacl) == [0, 0]
    @test hybridization(nacl) == [:none, :none]
end

@testset "rotatable" begin
    Phe = smilestomol("N[C@@H](CC1=CC=CC=C1)C(O)=O")
    @test rotatable_count(Phe) == 3
    KCl = smilestomol("[K+].[Cl-]")
    @test rotatable_count(KCl) == 0
    dipyridamole = smilestomol(
        "n3c(nc2c(nc(nc2N1CCCCC1)N(CCO)CCO)c3N4CCCCC4)N(CCO)CCO")
    @test rotatable_count(dipyridamole) == 12
    paclitaxel = smilestomol(
        "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
    )
    @test rotatable_count(paclitaxel) == 15
end

@testset "aromatic" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    phenol = smilestomol("C=1C=CC=CC=1O")  # not carbonyl O
    @test count(is_aromatic(phenol)) == 6

    furan = smilestomol("o1cccc1")
    @test count(is_aromatic(furan)) == 5

    quinone = smilestomol("C1(=O)C=CC(=O)C=C1")
    @test count(is_aromatic(quinone)) == 0

    tropone = smilestomol("C1(=O)C=CC=CC=C1")
    @test count(is_aromatic(tropone)) == 7

    azepine = smilestomol("N1C=CC=CC=C1")
    @test count(is_aromatic(azepine)) == 0

    thiopheneoxide = smilestomol("C1=CC=CS1=O")
    @test count(is_aromatic(thiopheneoxide)) == 0

    dihydroazaborine = smilestomol("C1=CC=CBN1")
    @test count(is_aromatic(dihydroazaborine)) == 6

    fulvene = smilestomol("C1=CC=CC1=C")
    @test count(is_aromatic(fulvene)) == 0

    cyclopropenyl = smilestomol("C1=C[C+]1")
    @test count(is_aromatic(cyclopropenyl)) == 3

    hexatriene = smilestomol("C=CC=CC=C")
    @test count(is_aromatic(hexatriene)) == 0

    borazine = smilestomol("B1NBNBN1")
    @test count(is_aromatic(borazine)) == 6

    methylenedioxybenzene = smilestomol("C1=CC=CC2=C1OCO2")
    @test count(is_aromatic(methylenedioxybenzene)) == 6

    naphthalene = smilestomol("C=1C=CC=C2C1C=CC=C2")
    @test count(is_aromatic(naphthalene)) == 10

    coumarin = smilestomol("C1=CC(=O)OC2=C1C=CC=C2")
    @test count(is_aromatic(coumarin)) == 10

    pyrene = smilestomol("C12=CC=C3C=CC=C4C=CC(C2=C34)=CC=C1")
    @test count(is_aromatic(pyrene)) == 16

    caffeine = smilestomol("CN1C=NC2=C1C(=O)N(C)C(=O)N2C")
    @test count(is_aromatic(caffeine)) == 9

    # Difficult aromaticities
    quinodimethane = smilestomol("C1=CC=CC(=C)C1=C")
    @test_broken count(is_aromatic(azulene)) == 10

    dihydronaphthalene = smilestomol("C=1CCC=C2C1C=CC=C2")
    @test_broken count(is_aromatic(dihydronaphthalene)) == 10

    azulene = smilestomol("C=1C=CC=2C1C=CC=CC2")
    @test_broken count(is_aromatic(azulene)) == 10

    pyromellitimide = smilestomol("C1(=O)NC(=O)C=2C1=CC=3C(=O)NC(=O)C=3C=2")
    @test count(is_aromatic(pyromellitimide)) == 6
    
    # TODO: how to deal with tautomerism
    isocyanurate = smilestomol("ClN1C(=O)N(Cl)C(=O)N(Cl)C1=O")
    # global_logger(default_logger)
end

@testset "properties" begin
    amide = smilestomol("CCC(=O)N")
    @test hydrogen_acceptor_count(amide) == 2
    @test hydrogen_donor_count(amide) == 1
    fluoro = smilestomol("CCN(CO)CF")
    @test hydrogen_acceptor_count(fluoro) == 3
    @test hydrogen_donor_count(fluoro) == 1
end

end # properties
