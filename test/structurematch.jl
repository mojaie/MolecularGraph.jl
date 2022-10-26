
using MolecularGraph:
    fastidentityfilter, fastsubstrfilter, emaptonmap


@testset "structurematch" begin

@testset "prefilter" begin
    BuOH = smilestomol("CCCO")
    butane = smilestomol("CCCC")
    propane = smilestomol("CCC")
    @test fastidentityfilter(butane, BuOH)
    @test !fastidentityfilter(propane, butane)
    @test fastsubstrfilter(butane, propane)
    @test !fastsubstrfilter(propane, butane)
end

@testset "structmatch" begin
    null = smilestomol("")
    @test !hasexactmatch(null, null)
    @test !hassubstructmatch(null, null)
    @test !hasnodeinducedmatch(null, null)
    @test !hasedgeinducedmatch(null, null)

    hexane = smilestomol("CCCCCC")
    iso = smilestomol("CC(C)CCC")
    iso2 = smilestomol("CCCC(C)C")
    cyclohexane = smilestomol("C1CCCCC1")
    @test !hasexactmatch(hexane, iso, prefilter=false)
    @test hasexactmatch(iso, iso2, prefilter=false)
    @test !hasexactmatch(cyclohexane, hexane, prefilter=false)
    @test hassubstructmatch(cyclohexane, hexane, prefilter=false)
    @test !hasnodeinducedmatch(cyclohexane, hexane, prefilter=false)
    @test hasedgeinducedmatch(cyclohexane, hexane, prefilter=false)

    tms = smilestomol("C[Si](C)(C)C")
    tsm = smilestomol("[Si]C([Si])([Si])[Si]")
    @test !hasexactmatch(tms, tsm, prefilter=false)
    @test !hasedgeinducedmatch(tms, tsm, prefilter=false)
    @test hasexactmatch(tms, tsm, prefilter=false, atommatcher=(m,q)->((m,q)->true))
    @test hasedgeinducedmatch(tms, tsm, prefilter=false, atommatcher=(m,q)->((m,q)->true))

    sulfide = smilestomol("CSSC")
    disconn = smilestomol("CS.SC")
    @test !hasexactmatch(sulfide, disconn, prefilter=false)
    @test hassubstructmatch(sulfide, disconn, prefilter=false)
    @test !hasnodeinducedmatch(sulfide, disconn, prefilter=false)
    @test hasedgeinducedmatch(sulfide, disconn, prefilter=false)

    # No edges
    nacl = smilestomol("[Na+].[Cl-]")
    na = smilestomol("[Na]")
    @test hassubstructmatch(nacl, na, prefilter=false)
    @test !hasedgeinducedmatch(nacl, na, prefilter=false)

    tetrahedrane = smilestomol("C12C3C1C23")
    fused = smilestomol("C1C2C1C2")
    spiro = smilestomol("C1CC12CC2")
    @test hassubstructmatch(tetrahedrane, fused, prefilter=false)
    @test !hassubstructmatch(tetrahedrane, spiro, prefilter=false)
    @test hasedgeinducedmatch(tetrahedrane, fused, prefilter=false)
    @test !hasedgeinducedmatch(tetrahedrane, spiro, prefilter=false)
end

@testset "connectedquery" begin
    hexane = smilestomol("CCCCCC")
    anyatom = smartstomol("*")
    @test hassubstructmatch(hexane, anyatom)

    aniline = smilestomol("C1=CC=CC=C1N")
    dieamine = smilestomol("CCNCC")
    priamine = smartstomol("[NX3;H2]")
    @test hassubstructmatch(aniline, priamine)
    @test !hassubstructmatch(dieamine, priamine)

    diether = smilestomol("CCOCC")
    phenol = smilestomol("c1ccccc1O")
    glycerol = smilestomol("OCC(O)CO")
    alcohol = smartstomol("[#6][OD]")
    @test !hassubstructmatch(diether, alcohol)
    @test hassubstructmatch(phenol, alcohol)
    @test hassubstructmatch(glycerol, alcohol)

    cyclopentane = smilestomol("C1CCCC1")
    pyrrole = smilestomol("N1C=CC=C1")
    aliphring = smartstomol("*@;!:*")
    @test hassubstructmatch(cyclopentane, aliphring)
    @test !hassubstructmatch(pyrrole, aliphring)

    triamine = smilestomol("CCN(CC)CC")
    acetamide = smilestomol("CC(=O)N")
    notamide = smartstomol(raw"[NX3;!$(NC=O)]")
    @test hassubstructmatch(triamine, notamide)
    @test !hassubstructmatch(acetamide, notamide)

    po1 = smilestomol("COOC")
    npo1 = smilestomol("COCOCOC")
    peroxide = smartstomol("[OX2][OX2]")
    @test hassubstructmatch(po1, peroxide)
    @test !hassubstructmatch(npo1, peroxide)

    pyridine = smilestomol("n1ccccc1")
    pyrrole = smilestomol("[nH]1cccc1")
    sixmem = smartstomol("[*r6]1[*r6][*r6][*r6][*r6][*r6]1")
    @test hassubstructmatch(pyridine, sixmem)
    @test !hassubstructmatch(pyrrole, sixmem)
end

@testset "disconnectedquery" begin
    hetero1 = smilestomol("CCN(O)[O-]")
    hetero2 = smilestomol("CCN=N")
    hetero3 = smilestomol("BrCCOC#N")
    hetero4 = smilestomol("CCS(=O)(=O)O")
    disconn = smartstomol("[#7,#8].[!#6].N")
    @test hassubstructmatch(hetero1, disconn)
    @test !hassubstructmatch(hetero2, disconn)
    @test hassubstructmatch(hetero3, disconn)
    @test !hassubstructmatch(hetero4, disconn)
end

@testset "querycontainment" begin
    tmse = smartstomol("O[Si](C)(C)C")
    tms = smartstomol("C[Si](C)C")
    @test hassubstructmatch(tmse, tms)
    @test !hassubstructmatch(tms, tmse)
    @test !hasexactmatch(tmse, tms)

    primary = smartstomol("[CH2][OH]")
    alcohol = smartstomol("C[OH]")
    @test hassubstructmatch(primary, alcohol)
    @test !hassubstructmatch(alcohol, primary)
    @test hasexactmatch(primary, alcohol)
    @test !hasexactmatch(alcohol, primary)
    @test !hasexactmatch(primary, alcohol, atommatcher=exactatommatch)

    pyrrole = smartstomol("[nH]1cccc1")
    pyridine = smartstomol("n1ccccc1")
    pyrrolidine = smartstomol("N1CCCC1")
    cyclicamine = smartstomol("[#7]1[#6][#6][#6][#6]1")
    @test !hassubstructmatch(pyridine, pyrrole)
    @test !hassubstructmatch(pyrrole, pyrrolidine)
    @test hassubstructmatch(pyrrole, cyclicamine)
    @test hassubstructmatch(pyrrolidine, cyclicamine)

    carbonylazide = smartstomol("O=CN=[N+]=[N-]")
    azide = smartstomol("N=[N+]=[N-]")
    hydrazine = smartstomol("N=N")
    stricthydrazine = smartstomol("[N+0]=[N+0]")
    @test hassubstructmatch(carbonylazide, azide)
    @test hassubstructmatch(azide, hydrazine)
    @test !hassubstructmatch(azide, stricthydrazine)

    halo = smartstomol("[#9,#17,#35]c1ccccc1")
    narrow = smartstomol("[#9,#17]c1ccccc1")
    broad = smartstomol("[#9,#17,#35,#53]c1ccccc1")
    broad2 = smartstomol("[#9,35#17,#35,#53]c1ccccc1")
    @test hasexactmatch(narrow, halo)
    @test !hasexactmatch(broad, halo)
    @test hasexactmatch(halo, broad)
    @test !hasexactmatch(halo, broad2)

    hetero = smartstomol("[!#6&!#7&!#8]1[#6][#6][#6][#6][#6]1")
    oxo = smartstomol("[#8]1[#6][#6][#6][#6][#6]1")
    sulfo = smartstomol("[#16]1[#6][#6][#6][#6][#6]1")
    thiopyrylium = smartstomol("[s+]1ccccc1")
    @test !hasexactmatch(oxo, hetero)
    @test hasexactmatch(sulfo, hetero)
    @test hasexactmatch(thiopyrylium, sulfo)
    @test hasexactmatch(thiopyrylium, hetero)

    diazocarbonyl1 = smartstomol(raw"[$(N=N=C~C=O)]")
    diazocarbonyl2 = smartstomol(raw"[$(N=N=C~C=O),$(N#N-C~C=O)]")
    @test hasexactmatch(diazocarbonyl1, diazocarbonyl2)
    @test !hasexactmatch(diazocarbonyl2, diazocarbonyl1)

    nested = smartstomol(raw"[$([CH]=[$(NOC)])]C=O")
    nestedor = smartstomol(raw"[$([CH]=[$(NOC),$(NO[Si])])]")
    @test hassubstructmatch(nested, nestedor)
    @test !hassubstructmatch(nestedor, nested)

    naryl = smartstomol(raw"OP(=O)(=[S,O])[$(Na)]")
    phos = smartstomol(raw"OP(=O)(=[S,O])N")
    phos2 = smartstomol(raw"[$([S,O]=PN)]")
    @test hasexactmatch(naryl, phos)
    @test !hasexactmatch(phos, naryl)
    @test hassubstructmatch(phos, phos2)

    tfas = smartstomol(raw"C(F)(F)(F)C(=O)S")
    tfmk1 = smartstomol(raw"[$(C(=O))]C(F)(F)(F)")
    tfmk2 = smartstomol(raw"[$(C(=O));!$(C-N);!$(C-O);!$(C-S)]C(F)(F)(F)")
    @test hassubstructmatch(tfas, tfmk1)
    @test !hassubstructmatch(tfas, tfmk2)

    quart = smartstomol(raw"[C+,Cl+,I+,P+,S+]")
    sulfonium = smartstomol(raw"[S+;X3;$(S-C);!$(S-[O;D1])]")
    cys = smartstomol(raw"NC(C=O)CS")
    @test hasexactmatch(sulfonium, quart)
    @test !hassubstructmatch(cys, sulfonium)
    @test !hassubstructmatch(cys, quart)

    halobenzene = smartstomol(raw"c1c([O;D1])c(-[Cl,Br,I])cc(-[Cl,Br,I])c1")
    halo = smartstomol(raw"c[F,Cl,Br,I]")
    @test hassubstructmatch(halobenzene, halo)

    pivoxil = smartstomol(raw"OCOC(=O)C([CH3])([CH3])[CH3]")
    ester = smartstomol(raw"[#6]-C(=O)O-[#6]")
    @test hassubstructmatch(pivoxil, ester)
end

@testset "node matching" begin
    function nodematch(mol, query, idx=1, len=1)
        afunc = atommatch(mol, query)
        bfunc = bondmatch(mol, query)
        matches = edgesubgraph_isomorphisms(mol, query, nodematcher=afunc, edgematcher=bfunc)
        emaps = collect(matches)
        @test length(emaps) == len
        emap = emaps[idx]
        return emaptonmap(emap, mol, query)
    end
    function single(v)   # like `only` on Julia 1.4+
        length(v) == 1 || error("expected a single entry")
        return v[1]
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
    @test single.(nmap) == [8,7,6,5,10,9,1,2,3,4]  # indices in mol corresponding to atoms in query
    # 1-Acenaphthenone, PubChem CID75229
    mol = smilestomol("C1C2=CC=CC3=C2C(=CC=C3)C1=O")
    nmap = nodematch(mol, query)
    @test single.(nmap) == [11,10,9,8,7,6,2,1,12,13]
    # Ninhydrin, PubChem CID10236 (2 matches)
    mol = smilestomol("C1=CC=C2C(=C1)C(=O)C(C2=O)(O)O")
    nmaps = ([1,2,3,4,5,6,7,9,10,11], [2,1,6,5,4,3,10,9,7,8])
    nmap = nodematch(mol, query, 1, 2)
    snmap = single.(nmap)
    @test snmap ∈ nmaps
    nmap = nodematch(mol, query, 2, 2)
    snmap = single.(nmap)
    @test snmap ∈ nmaps

    # Ambiguous match
    query = smartstomol("[H][H]")
    mol = smilestomol("[H][H]")
    nmap = nodematch(mol, query)
    @test nmap == [[1,2], [1,2]]
end

@testset "MCS" begin
    null = smilestomol("")
    @test size(disconnectedmcis(null, null)) == 0

    BuOH = smilestomol("CCCO")
    butane = smilestomol("CCCC")
    @test size(disconnectedmcis(butane, BuOH)) == 3
    @test size(disconnectedmces(butane, BuOH)) == 2

    hexane = smilestomol("CCCCCC")
    cyclohexane = smilestomol("C1CCCCC1")
    @test size(disconnectedmcis(hexane, cyclohexane)) == 5
    @test size(disconnectedmces(hexane, cyclohexane)) == 5

    tms = smilestomol("C[Si](C)(C)C")
    tsm = smilestomol("[Si]C([Si])([Si])[Si]")
    @test size(disconnectedmcis(tms, tsm)) == 2
    @test size(disconnectedmces(tms, tsm)) == 1

    sulfide = smilestomol("CSSC")
    disconn = smilestomol("CS.SC")
    @test size(connectedmcis(sulfide, disconn)) == 2
    @test size(connectedmces(sulfide, disconn)) == 1
    @test size(disconnectedmcis(sulfide, disconn)) == 3
    @test size(disconnectedmces(sulfide, disconn)) == 2

    # Delta-Y
    cyclopropane = smilestomol("C1CC1")
    isopropane = smilestomol("CC(C)C")
    @test size(disconnectedmcis(cyclopropane, isopropane)) == 2
    @test size(disconnectedmces(cyclopropane, isopropane)) == 2

    tetrahedrane = smilestomol("C12C3C1C23")
    fused = smilestomol("C1C2C1C2")
    spiro = smilestomol("C1CC12CC2")
    @test size(disconnectedmcis(tetrahedrane, fused)) == 3
    @test size(disconnectedmcis(tetrahedrane, spiro)) == 3
    @test size(disconnectedmces(tetrahedrane, fused)) == 5
    @test size(disconnectedmces(tetrahedrane, spiro)) == 4

    # Topological constraint
    sulfide = smilestomol("CSSC")
    disconn = smilestomol("CS.SC")
    @test size(tcmcis(sulfide, disconn)) == 2
    @test size(tcmces(sulfide, disconn)) == 1

    cid6437877 = smilestomol("CC1=C(SC=N1)C=CC2=C(N3C(C(C3=O)NC(=O)C(=NOC)C4=CSC(=N4)N)SC2)C(=O)OCOC(=O)C(C)(C)C")  # Cefditoren pivoxil
    cid5481173 = smilestomol("CC(C)(C(=O)O)ON=C(C1=CSC(=N1)N)C(=O)NC2C3N(C2=O)C(=C(CS3)C[N+]4=CC=CC=C4)C(=O)[O-]")  # Ceftazidime
    @test size(tcmcis(cid6437877, cid5481173, tolerance=1)) == 27
    @test size(tcmces(cid6437877, cid5481173, tolerance=1)) == 28
    @test size(tcmcis(cid6437877, cid5481173, diameter=8)) == 18
    @test size(tcmces(cid6437877, cid5481173, diameter=8)) == 20
end

end # structurematch
