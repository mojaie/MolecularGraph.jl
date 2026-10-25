
@testset "properties" begin

@testset "descriptor" begin
    desc = MolDescriptor{Int}()
    push!(desc.sssr, collect(2:7), collect(9:14))
    remap!(
        Val(:sssr), desc, [1, 2, 14, 4, 5, 13, 7, 8, 9, 10, 11, 12],
        Edge.([(1, 2), (2, 3)]) # Not used
    )
    @test length(desc.sssr) == 1
    @test desc.sssr[1] == [9, 10, 11, 12, 6, 3]
    dump = JSON.json(desc.sssr)
    @test JSON.parse(dump, Vector{Vector{Int}}) == desc.sssr
    cp = copy(desc)
    desc.sssr[1][3] == 13
    cp.sssr[1][3] == 11
end

@testset "topology" begin
    cubane = smilestomol("C12C3C4C1C1C4C3C12")
    @test all(is_in_ring(cubane))
    @test smallest_ring(cubane)[4] == 4
    biphenyl = smilestomol("C1CCCCC1C1CCCCC1")
    @test ring_count(biphenyl)[12] == 1
    @test is_in_ring(biphenyl)[6]
    @test is_in_ring(biphenyl)[7]
    @test !is_edge_in_ring(biphenyl)[edge_rank(edge_rank(biphenyl), 6, 7)]
    rem_vertex!(biphenyl, 5)
    # TODO: better cycle equivalence check
    @test sum(only(sssr(biphenyl))) == 50  # remapped, [7, 8, 9, 10, 11, 5]

    subpyr = smartstomol("c1[cX3]nccc1")
    @test length(sssr(subpyr)) == 1
end

@testset "elements" begin
    sodiumsulfate = smilestomol("O=S(=O)([O-])[O-].[Na+].[Na+]")
    @test atom_symbol(sodiumsulfate) == [:O, :S, :O, :O, :O, :Na, :Na]
    @test bond_order(sodiumsulfate) == [2, 2, 1, 1]
    @test atom_charge(sodiumsulfate) == [0, 0, 0, -1, -1, 1, 1]
    @test degree(sodiumsulfate) == [1, 4, 1, 1, 1, 0, 0]

    isocyanurate = smilestomol("ClN1C(=O)N(Cl)C(=O)N(Cl)C1=O")
    @test atom_symbol(isocyanurate) == [
        :Cl, :N, :C, :O, :N, :Cl, :C, :O, :N, :Cl, :C, :O]
    @test bond_order(isocyanurate) == [
        1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2]
    @test all(atom_charge(isocyanurate) .== 0)
    @test degree(isocyanurate) == [1, 3, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1]

    aldehyde = smartstomol("[CH]=O")
    @test_throws MethodError atom_symbol(aldehyde)
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
    @test molecular_formula(smilestomol("[Na+].[Cl-]")) == "ClNa"
    @test molecular_formula(smilestomol("ClC(Cl)(Cl)Cl")) == "CCl4"
end

@testset "valence" begin
    atoms = smilestomol("B.C.N.O.F.[SiH4].P.S.Cl.[AsH3].[SeH2].Br.I")
    heavy = findall(!=(:H), atom_symbol(atoms))
    @test valence(atoms)[heavy] == [3, 4, 3, 2, 1, 4, 3, 2, 1, 3, 2, 1, 1]
    @test lone_pair(atoms)[heavy] == [0, 0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3, 3]
    @test total_hydrogens(atoms)[heavy] == [3, 4, 3, 2, 1, 4, 3, 2, 1, 3, 2, 1, 1]

    # bracket atoms have no implicit hydrogens
    bare = smilestomol("[Si].[As].[Se].[C].[Cl]")
    @test valence(bare) == [0, 0, 0, 0, 0]
    @test implicit_hydrogens(bare) == [0, 0, 0, 0, 0]
    @test atom_valence(bare) == [0, 0, 0, 0, 0]
    radical = smilestomol("[CH3]")
    @test atom_valence(radical)[1] == 3
    @test molecular_formula(radical) == "CH3"
    noval = smilestomol("[C@@H](F)(Cl)Br.[NH4+].[O-]C.c1cc[nH]c1")
    @test all(==(-1), atom_valence(noval))  # not set if the bond order is enough
    # non-organic elements always have the specified valence
    metal = smilestomol("[Sn].C[Sn](C)C.[Na+].[SnH4]")
    @test atom_valence(metal)[findall(!=(:H), atom_symbol(metal))] == [0, -1, 3, -1, -1, 0, 4]
    @test atom_valence(smilestomol("[2H]C.[H][H]")) == [-1, -1, -1, -1]
    # aromatic bracket atoms
    phenyl = smilestomol("c1cc[c]cc1")
    @test atom_valence(phenyl) == [-1, -1, -1, 3, -1, -1]
    @test phenyl[4].isaromatic
    @test implicit_hydrogens(phenyl)[4] == 0
    @test molecular_formula(phenyl) == "C6H5"
    # apparent valence of aromatic atoms before kekulization
    nokek(m) = (MolecularGraph.default_atom_charge!(m); MolecularGraph.default_bond_order!(m);
        apparent_valence!(m); valence!(m))
    indole = smilestomol("c1ccc2[nH]ccc2c1"; on_update=nokek)
    @test apparent_valence(indole) == [3, 3, 3, 4, 3, 1, 3, 3, 4, 3]
    @test implicit_hydrogens(indole) == [1, 1, 1, 0, 0, 0, 1, 1, 0, 1]
    @test apparent_valence(smilestomol("n1cccc1"; on_update=nokek)) == [3, 3, 3, 3, 3]
    @test apparent_valence(smilestomol("O=c1cccc[nH]1"; on_update=nokek)) == [2, 4, 3, 3, 3, 3, 3, 1]
    @test apparent_valence(smilestomol("c1ccn(=O)cc1"; on_update=nokek))[4] == 5
    # not changed after kekulization
    @test apparent_valence(smilestomol("c1ccc2[nH]ccc2c1")) == apparent_valence(indole)
    aromnoval = smilestomol("c1cc[cH]cc1.c1cc[nH]c1.c1cc[n]cc1.c1cc[nH+]cc1.c1cc[o]c1.C[n]1cccc1.O=[c]1cccc[nH]1")
    @test all(==(-1), atom_valence(aromnoval))
    @test molecular_formula(aromnoval) == "C34H38N5O2"
    @test occursin(r"Sn  0  0  0  0  0 15", printv2mol(smilestomol("[Sn]")))

    # carbene
    carbene = smilestomol("[CH2]")
    @test valence(carbene) == [2, 1, 1]
    @test implicit_hydrogens(carbene) == [0, 0, 0]
    @test lone_pair(carbene)[1] == 1
    @test hybridization(carbene)[1] === :sp2
    @test molecular_formula(carbene) == "CH2"
    triplet = smilestomol("C")
    triplet[1] = SMILESAtom(symbol=:C, multiplicity=3, valence=2)
    @test molecular_formula(triplet) == "CH2"
    @test lone_pair(triplet) == [0]
    dichlorocarbene = smilestomol("Cl[C]Cl")
    @test valence(dichlorocarbene) == [1, 2, 1]
    @test lone_pair(dichlorocarbene) == [3, 1, 3]

    # pentavalent nitrogen and nitro group
    nitro = smilestomol("CN(=O)=O")
    @test valence(nitro) == [4, 5, 2, 2]
    @test lone_pair(nitro) == [0, 0, 2, 2]
    @test hybridization(nitro)[2] === :sp2
    pentan = smilestomol("CN")
    pentan[2] = SMILESAtom(symbol=:N, valence=5)
    @test valence(pentan) == [4, 5]
    @test implicit_hydrogens(pentan) == [3, 4]
    @test lone_pair(pentan) == [0, 0]
    @test molecular_formula(pentan) == "CH7N"
    # specified valence lower than the bond order
    low = smilestomol("CC(C)C")
    low[2] = SMILESAtom(symbol=:C, valence=2)
    @test valence(low)[2] == 3

    pyridineoxide = smilestomol("[N+]1([O-])=CC=CC=C1")
    @test valence(pyridineoxide) == [4, 1, 4, 4, 4, 4, 4]
    @test lone_pair(pyridineoxide) == [0, 3, 0, 0, 0, 0, 0]
    @test implicit_hydrogens(pyridineoxide) == [0, 0, 1, 1, 1, 1, 1]

    trifluoroborate = smilestomol("C[B-](F)(F)F")
    @test valence(trifluoroborate) == [4, 4, 1, 1, 1]
    @test lone_pair(trifluoroborate) == [0, 0, 3, 3, 3]
    @test implicit_hydrogens(trifluoroborate) == [3, 0, 0, 0, 0]

    etmacl = smilestomol("[H]C([H])([H])C[Mg][Cl]")
    @test valence(etmacl) == [1, 4, 1, 1, 4, 2, 1]
    @test lone_pair(etmacl) == [0, 0, 0, 0, 0, 0, 3]
    @test implicit_hydrogens(etmacl) == [0, 0, 0, 0, 2, 0, 0]
    @test explicit_hydrogens(etmacl) == [0, 3, 0, 0, 0, 0, 0]
    @test total_hydrogens(etmacl) == [0, 3, 0, 0, 2, 0, 0]
    @test connectivity(etmacl) == [1, 4, 1, 1, 4, 2, 1]
end

@testset "hybridization" begin
    alkyne = smilestomol("CC=CC#C")
    @test pi_electron(alkyne) == [0, 1, 1, 2, 2]
    @test hybridization(alkyne) == [:sp3, :sp2, :sp2, :sp, :sp]

    nitrile = smilestomol("C#N")
    @test pi_electron(nitrile) == [2, 2]
    @test hybridization(nitrile) == [:sp, :sp]

    carboxylate = smilestomol("CC(=O)[O-]")
    @test pi_electron(carboxylate) == [0, 1, 1, 2]
    @test hybridization(carboxylate) == [:sp3, :sp2, :sp2, :sp2]

    ester = smilestomol("CC(=O)OC")
    @test pi_electron(ester) == [0, 1, 1, 2, 0]
    @test hybridization(ester) == [:sp3, :sp2, :sp2, :sp2, :sp3]

    anilinium = smilestomol("C1=CC=CC=C1[NH3+]")  # bracket atoms have no implicit hydrogens
    @test pi_electron(anilinium)[1:7] == [1, 1, 1, 1, 1, 1, 0]
    @test hybridization(anilinium)[1:7] == [:sp2, :sp2, :sp2, :sp2, :sp2, :sp2, :sp3]

    quat = smilestomol("C[N+](C)(C)C(=O)N")
    @test pi_electron(quat) == [0, 0, 0, 0, 1, 1, 2]
    @test hybridization(quat) == [:sp3, :sp3, :sp3, :sp3, :sp2, :sp2, :sp2]

    azide = smilestomol("CC(=O)N=[N+]=[N-]")
    @test pi_electron(azide) == [0, 1, 1, 1, 2, 1]
    @test hybridization(azide) == [:sp3, :sp2, :sp2, :sp2, :sp, :sp2]

    borane = smilestomol("B")
    @test pi_electron(borane) == [0]
    @test hybridization(borane) == [:sp2]

    carbocation = smilestomol("[SiH3][C+]([SiH3])[SiH3]")
    heavy = findall(!=(:H), atom_symbol(carbocation))
    @test pi_electron(carbocation)[heavy] == [0, 0, 0, 0]
    @test hybridization(carbocation)[heavy] == [:sp3, :sp2, :sp3, :sp3]

    pyrrole = smilestomol("C1=CC=CN1")
    @test pi_electron(pyrrole) == [1, 1, 1, 1, 2]
    @test hybridization(pyrrole) == [:sp2, :sp2, :sp2, :sp2, :sp2]

    thiophene = smilestomol("S1C=CC=C1")
    @test pi_electron(thiophene) == [2, 1, 1, 1, 1]
    @test hybridization(thiophene) == [:sp2, :sp2, :sp2, :sp2, :sp2]

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

    thiophene = smilestomol("S1C=CC=C1")
    @test count(is_aromatic(thiophene)) == 5

    quinone = smilestomol("C1(=O)C=CC(=O)C=C1")
    @test count(is_aromatic(quinone)) == 0

    tropone = smilestomol("C1(=O)C=CC=CC=C1")
    @test count(is_aromatic(tropone)) == 7

    azepine = smilestomol("N1C=CC=CC=C1")
    @test count(is_aromatic(azepine)) == 0

    borepin = smilestomol("B1C=CC=CC=C1")
    @test count(is_aromatic(borepin)) == 7

    thiopheneoxide = smilestomol("C1=CC=CS1=O")
    @test count(is_aromatic(thiopheneoxide)) == 0

    dihydroazaborine = smilestomol("C1=CC=CBN1")
    @test count(is_aromatic(dihydroazaborine)) == 6

    fulvene = smilestomol("C1=CC=CC1=C")
    @test count(is_aromatic(fulvene)) == 0

    cyclopropenyl = smilestomol("C1=C[CH+]1")
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

    imidazo12a = smilestomol("C1=CN=C2N1C=CC=C2")  # imidazo[1,2-a]pyridine
    @test count(is_aromatic(imidazo12a)) == 9

    coe = smilestomol(
        "C1=CC=C(C=C1)CC2=C3N=C(C(=O)N3C=C(N2)C4=CC=C(C=C4)O)CC5=CC=C(C=C5)O")  # coelenterazine
    @test count(is_aromatic(coe)) == 27

    quinodimethane = smilestomol("C1=CC=CC(=C)C1=C")
    @test count(is_aromatic(quinodimethane)) == 0

    azulene = smilestomol("C=1C=CC=2C1C=CC=CC2")
    @test count(is_aromatic(azulene)) == 10

    pyromellitimide = smilestomol("C1(=O)NC(=O)C=2C1=CC=3C(=O)NC(=O)C=3C=2")
    @test count(is_aromatic(pyromellitimide)) == 6

    dihydroanth = smilestomol("C2=CC=CC3=CC1=CCCC=C1C=C23")  # 2,3-dihydroanthracene 
    @test count(is_aromatic(dihydroanth)) == 0

    cyclononen = smilestomol("C14=CC=CC1=CC2=CC=CC2=CC3=CC=CC3=C4")
    @test count(is_aromatic(cyclononen)) == 18

    indacene = smilestomol("C1=CC2=CC3=CC=CC3=CC2=C1")  # s-indacene
    @test count(is_aromatic(indacene)) == 0

    isocyanurate = smilestomol("ClN1C(=O)N(Cl)C(=O)N(Cl)C1=O")
    @test count(is_aromatic(isocyanurate)) == 6

    c20 = smilestomol(
        "C=12C3=C4C5=C6C7=C(C=15)C1=C2C2=C5C1=C7C1=C5C(=C32)C4=C61")  # C20 fullerene
    @test count(is_aromatic(c20)) == 0

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

@testset "update" begin
    pr = smilestomol("CCC")
    add_vertex!(pr, SMILESAtom(:O))
    add_edge!(pr, 3, 4, SMILESBond(2))
    @test sum(bond_order(pr)) == 4
    add_vertex!(pr, SMILESAtom(:N, 1))
    add_edge!(pr, 1, 5, SMILESBond())
    @test sum(atom_charge(pr)) == 1
    rem_vertex!(pr, 4)
    @test sum(bond_order(pr)) == 3
    pr[4] = SMILESAtom(:N)
    @test lone_pair(pr)[4] == 1
    rem_vertices!(pr, [2, 3])
    @test sum(valence(pr)) == 7
end

@testset "scaffold_fragments" begin
    scfs, dag = MolecularGraph.scaffold_fragments(smilestomol("CCO").graph)  # no scaffold
    @test isempty(scfs) && nv(dag) == 0
    scfs, dag = MolecularGraph.scaffold_fragments(smilestomol("c1ccccc1Cc1ccncc1").graph)
    @test length(scfs) == nv(dag) >= 1
end

end # properties
