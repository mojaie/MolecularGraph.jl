
using MolecularGraph:
    exact_match_prefilter, substruct_match_prefilter


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

@testset "querycontainment" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))

    # no properties
    fmr = smartstomol("*1***1")
    three = smartstomol("***")
    @test has_substruct_match(fmr, three)
    @test !has_substruct_match(three, fmr)

    tmse = smartstomol("O[Si](C)(C)C")
    tms = smartstomol("C[Si](C)C")
    @test has_substruct_match(tmse, tms)
    @test !has_exact_match(tmse, tms)

    primary = smartstomol("[CH2][OH]")
    alcohol = smartstomol("C[OH]")
    @test has_exact_match(primary, alcohol)
    @test !has_exact_match(alcohol, primary)

    pyrrole = smartstomol("[nH]1cccc1")
    pyridine = smartstomol("n1ccccc1")
    pyrrolidine = smartstomol("N1CCCC1")
    cyclicamine = smartstomol("[#7]1[#6][#6][#6][#6]1")
    @test !has_substruct_match(pyridine, pyrrole)
    @test !has_substruct_match(pyrrole, pyrrolidine)
    @test has_exact_match(pyrrole, cyclicamine)
    @test has_exact_match(pyrrolidine, cyclicamine)

    carbonylazide = smartstomol("O=CN=[N+]=[N-]")
    azide = smartstomol("N=[N+]=[N-]")
    hydrazine = smartstomol("N=N")
    stricthydrazine = smartstomol("[N+0]=[N+0]")
    @test has_substruct_match(carbonylazide, azide)
    @test has_substruct_match(azide, hydrazine)
    @test !has_substruct_match(azide, stricthydrazine)

    halo = smartstomol("[#9,#17,#35]c1ccccc1")
    narrow = smartstomol("[#9,#17]c1ccccc1")
    broad = smartstomol("[#9,#17,#35,#53]c1ccccc1")
    broad2 = smartstomol("[#9,35#17,#35,#53]c1ccccc1")
    @test has_exact_match(narrow, halo)
    @test !has_exact_match(broad, halo)
    @test has_exact_match(halo, broad)
    @test !has_exact_match(halo, broad2)

    hetero = smartstomol("[!#6&!#7&!#8]1[#6][#6][#6][#6][#6]1")
    oxo = smartstomol("[#8]1[#6][#6][#6][#6][#6]1")
    sulfo = smartstomol("[#16]1[#6][#6][#6][#6][#6]1")
    thiopyrylium = smartstomol("[s+]1ccccc1")
    @test !has_exact_match(oxo, hetero)
    @test has_exact_match(sulfo, hetero)
    @test has_exact_match(thiopyrylium, sulfo)
    @test has_exact_match(thiopyrylium, hetero)

    diazocarbonyl1 = smartstomol(raw"[$(N=N=C~C=O)]")
    diazocarbonyl2 = smartstomol(raw"[$(N=N=C~C=O),$(N#N-C~C=O)]")
    @test has_exact_match(diazocarbonyl1, diazocarbonyl2)
    @test !has_exact_match(diazocarbonyl2, diazocarbonyl1)

    nested = smartstomol(raw"[$([CH]=[$(NOC)])]C=O")
    nestedor = smartstomol(raw"[$([CH]=[$(NOC),$(NO[Si])])]")
    @test has_substruct_match(nested, nestedor)
    @test !has_substruct_match(nestedor, nested)

    naryl = smartstomol(raw"OP(=O)(=[S,O])[$(Na)]")
    phos = smartstomol(raw"OP(=O)(=[S,O])N")
    phos2 = smartstomol(raw"[$([S,O]=PN)]")
    @test has_exact_match(naryl, phos)
    @test !has_exact_match(phos, naryl)
    @test !has_substruct_match(phos, phos2)

    tfas = smartstomol(raw"C(F)(F)(F)C(=O)S")
    tfmk1 = smartstomol(raw"PS[$(C(=O))](=O)C(F)(F)(F)")
    tfmk2 = smartstomol(raw"[$(C(=O));!$(C-N);!$(C-O);!$(C-S)]C(F)(F)(F)")
    @test has_substruct_match(tfmk1, tfas)
    @test !has_substruct_match(tfas, tfmk2)

    quart = smartstomol(raw"[C+,Cl+,I+,P+,S+]")
    sulfonium = smartstomol(raw"[S+;X3;$(S-C);!$(S-[O;D1])]")
    cys = smartstomol(raw"NC(C=O)CS")
    @test has_exact_match(sulfonium, quart)
    @test !has_substruct_match(cys, sulfonium)
    @test !has_substruct_match(cys, quart)

    halobenzene = smartstomol(raw"c1c([O;D1])c(-[Cl,Br,I])cc(-[Cl,Br,I])c1")
    halo = smartstomol(raw"c[F,Cl,Br,I]")
    @test has_substruct_match(halobenzene, halo)

    pivoxil = smartstomol(raw"OCOC(=O)C([CH3])([CH3])[CH3]")
    ester = smartstomol(raw"[#6]C(=O)O[#6]")
    ester2 = smartstomol(raw"[#6]-C(=O)O-[#6]")  # MLSMR structural alert Ester
    @test has_substruct_match(pivoxil, ester)
    @test !has_substruct_match(pivoxil, ester2)

    carbon = smartstomol("C")
    # BMS structural alert 'contains metal'
    # Just for int overflow check (set max property size cutoff, do not generate 2^61 truthtable)
    metals = smartstomol(raw"[$([Ru]),$([Rh]),$([Se]),$([Pd]),$([Sc]),$([Bi]),$([Sb]),$([Ag]),$([Ti]),$([Al]),$([Cd]),$([V]),$([In]),$([Cr]),$([Sn]),$([Mn]),$([La]),$([Fe]),$([Er]),$([Tm]),$([Yb]),$([Lu]),$([Hf]),$([Ta]),$([W]),$([Re]),$([Co]),$([Os]),$([Ni]),$([Ir]),$([Cu]),$([Zn]),$([Ga]),$([Ge]),$([As]),$([Y]),$([Zr]),$([Nb]),$([Ce]),$([Pr]),$([Nd]),$([Sm]),$([Eu]),$([Gd]),$([Tb]),$([Dy]),$([Ho]),$([Pt]),$([Au]),$([Hg]),$([Tl]),$([Pb]),$([Ac]),$([Th]),$([Pa]),$([Mo]),$([U]),$([Tc]),$([Te]),$([Po]),$([At])]")
    @test !has_substruct_match(carbon, metals)

    # global_logger(default_logger)
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
