
@testset "structurematch" begin

@testset "prefilter" begin
    BuOH = smilestomol("CCCO")
    butane = smilestomol("CCCC")
    propane = smilestomol("CCC")
    @test exact_match_prefilter(butane, BuOH)
    @test !exact_match_prefilter(propane, butane)
    @test substruct_match_prefilter(butane, propane)
    @test !substruct_match_prefilter(propane, butane)
end

@testset "attribute_match" begin
    mol = smilestomol("[N+][14C@@H](C)C=O")
    @test string(atom_symbol(mol)[6]) == "O"  # qeq(:symbol, "O")
    @test string(~is_aromatic(mol)[2]) == "true"  # qtrue(:isaromatic)
    @test string(atom_charge(mol)[1]) == "1"  # qeq(:charge, "1" )
    @test string([atom_mass(props(mol, i)) for i in vertices(mol)][2]
        ) == "14"  # qeq(:mass => "14")
    @test string([atom_mass(props(mol, i)) for i in vertices(mol)][4]
        ) == "nothing"  # qeq(:mass => "nothing")
    @test string(connectivity(mol)[2]) == "4"  # qeq(:connectivity, "4")
    @test string(degree(mol)[2]) == "4"  # qeq(:degree, "4")
    @test string(valence(mol)[6]) == "2"  # qeq(:valence, "2")
    @test string(total_hydrogens(mol)[2]) == "1"  # qeq(:total_hydrogens, "1")
    @test string(smallest_ring(mol)[3]) == "0"  # qeq(:smallest_ring, "0")
    @test string(ring_count(mol)[4]) == "0"  # qeq(:ring_count, "0")

    @test string(bond_order(mol)[5]) == "2"  # qeq(:order, "2")
    @test string(~is_edge_in_ring(mol)[4]) == "true"  # qnot(), qeq(:is_in_ring, "true")
    @test string(~is_edge_aromatic(mol)[4]) == "true"  # qnot(), qeq(:isaromatic, "true")
end

@testset "structmatch" begin
    null = smilestomol("")
    @test !has_exact_match(null, null)
    @test !has_substruct_match(null, null)

    iso = smilestomol("CC(C)CCC")
    iso2 = smilestomol("CCCC(C)C")
    hexane = smilestomol("CCCCCC")
    cyclohexane = smilestomol("C1CCCCC1")
    @test has_exact_match(iso, iso2)
    @test !has_exact_match(hexane, iso)
    @test has_substruct_match(cyclohexane, hexane)

    tms = smilestomol("C[Si](C)(C)C")
    tsm = smilestomol("[Si]C([Si])([Si])[Si]")
    @test !has_exact_match(tms, tsm)

    sulfide = smilestomol("CSSC")
    disconn = smilestomol("CS.SC")
    @test !has_exact_match(sulfide, disconn)
    @test has_substruct_match(sulfide, disconn)

    nacl = smilestomol("[Na+].[Cl-]")  # No edges
    na = smilestomol("[Na]")
    @test has_substruct_match(nacl, na)

    tetrahedrane = smilestomol("C12C3C1C23")
    fused = smilestomol("C1C2C1C2")
    spiro = smilestomol("C1CC12CC2")
    @test has_substruct_match(tetrahedrane, fused)
    @test !has_substruct_match(tetrahedrane, spiro)

    # geometry

    cumene = smilestomol("c1ccccc1C(C)C")
    isobutane = smilestomol("CC(C)C")
    propane = smilestomol("CCC")
    @test !has_substruct_match(cumene, isobutane)
    @test has_substruct_match(cumene, propane)

    thiazole = smilestomol("s1cncc1N")
    ethylene = smilestomol("C=CN")
    @test has_substruct_match(thiazole, ethylene)

    diazo1 = smilestomol("C=[N+]=[N-]")  # sp2, sp, sp2
    diazo2 = smilestomol("[C-][N+]#N")  # sp3, sp, sp
	nitrogen = smilestomol("N#N")  # sp, sp
    @test !has_substruct_match(diazo1, nitrogen)
    @test has_substruct_match(diazo2, nitrogen)

    sulfone = smilestomol("CS(=O)(=O)N")
	sulfoxide = smilestomol("CS(=O)N")
    @test has_substruct_match(sulfone, sulfoxide)
end

@testset "connectedquery" begin
    hexane = smilestomol("CCCCCC")
    anyatom = smartstomol("*")
    @test has_substruct_match(hexane, anyatom)

    aniline = smilestomol("C1=CC=CC=C1N")
    dieamine = smilestomol("CCNCC")
    priamine = smartstomol("[NX3;H2]")
    @test has_substruct_match(aniline, priamine)
    @test !has_substruct_match(dieamine, priamine)

    diether = smilestomol("CCOCC")
    phenol = smilestomol("c1ccccc1O")
    glycerol = smilestomol("OCC(O)CO")
    alcohol = smartstomol("[#6][OD]")
    @test !has_substruct_match(diether, alcohol)
    @test has_substruct_match(phenol, alcohol)
    @test has_substruct_match(glycerol, alcohol)

    cyclopentane = smilestomol("C1CCCC1")
    pyrrole = smilestomol("N1C=CC=C1")
    aliphring = smartstomol("*@;!:*")
    @test has_substruct_match(cyclopentane, aliphring)
    @test !has_substruct_match(pyrrole, aliphring)

    triamine = smilestomol("CCN(CC)CC")
    acetamide = smilestomol("CC(=O)N")
    notamide = smartstomol(raw"[NX3;!$(NC=O)]")
    @test has_substruct_match(triamine, notamide)
    @test !has_substruct_match(acetamide, notamide)

    po1 = smilestomol("COOC")
    npo1 = smilestomol("COCOCOC")
    peroxide = smartstomol("[OX2][OX2]")
    @test has_substruct_match(po1, peroxide)
    @test !has_substruct_match(npo1, peroxide)

    pyridine = smilestomol("n1ccccc1")
    pyrrole = smilestomol("[nH]1cccc1")
    sixmem = smartstomol("[*r6]1[*r6][*r6][*r6][*r6][*r6]1")
    @test has_substruct_match(pyridine, sixmem)
    @test !has_substruct_match(pyrrole, sixmem)
end

@testset "disconnectedquery" begin
    hetero1 = smilestomol("CCN(O)[O-]")
    hetero2 = smilestomol("CCN=N")
    hetero3 = smilestomol("BrCCOC#N")
    hetero4 = smilestomol("CCS(=O)(=O)O")
    disconn = smartstomol("[#7,#8].[!#6].N")
    @test has_substruct_match(hetero1, disconn)
    @test !has_substruct_match(hetero2, disconn)
    @test has_substruct_match(hetero3, disconn)
    @test !has_substruct_match(hetero4, disconn)
end

@testset "node matching" begin
    function nodematch(mol, query, idx=1, len=1)
        afunc = vmatchgen(mol, query)
        bfunc = ematchgen(mol, query)
        matches = edgesubgraph_isomorphisms(mol.graph, query.graph, vmatch=afunc, ematch=bfunc)
        emaps = collect(matches)
        @test length(emaps) == len
        emap = emaps[idx]
        return emaptonmap(emap, mol, query)
    end
    # 6-membered carbon ring fused to a 5-membered carbon ring with an attached oxygen
    query = smartstomol("[#6]~1~[#6]~[#6]~[#6]~2~[#6](~[#6]~1)~[#6]~[#6]~[#6]~2~[#8]")
    # The carbon numbering looks like this:
    #         3     9
    #      2     4     8
    #      1     5     7
    #         6
    # 1-Indanone, PubChem CID6735
    mol = smilestomol("C1CC(=O)C2=CC=CC=C21")
    nmap = nodematch(mol, query)
    @test only.(nmap) == [8,7,6,5,10,9,1,2,3,4]  # indices in mol corresponding to atoms in query
    # 1-Acenaphthenone, PubChem CID75229
    mol = smilestomol("C1C2=CC=CC3=C2C(=CC=C3)C1=O")
    nmap = nodematch(mol, query)
    @test only.(nmap) == [11,10,9,8,7,6,2,1,12,13]
    # Ninhydrin, PubChem CID10236 (2 matches)
    mol = smilestomol("C1=CC=C2C(=C1)C(=O)C(C2=O)(O)O")
    nmaps = ([1,2,3,4,5,6,7,9,10,11], [2,1,6,5,4,3,10,9,7,8])
    nmap = nodematch(mol, query, 1, 2)
    snmap = only.(nmap)
    @test snmap ∈ nmaps
    nmap = nodematch(mol, query, 2, 2)
    snmap = only.(nmap)
    @test snmap ∈ nmaps

    # Ambiguous match
    query = smartstomol("[H][H]")
    mol = smilestomol("[H][H]")
    nmap = nodematch(mol, query)
    @test nmap == [[1,2], [1,2]]
end

end # structurematch
