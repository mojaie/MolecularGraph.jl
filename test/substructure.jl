
using MolecularGraph:
    fastidentityfilter, fastsubstrfilter,
    atommatch, bondmatch, emaptonmap


@testset "structmatch" begin

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
    @test !isstructmatch(null, null, :exact, prefilter=false)
    @test !isstructmatch(null, null, :edgeinduced, prefilter=false)

    hexane = smilestomol("CCCCCC")
    iso = smilestomol("CC(C)CCC")
    iso2 = smilestomol("CCCC(C)C")
    cyclohexane = smilestomol("C1CCCCC1")
    @test !isstructmatch(hexane, iso, :exact, prefilter=false)
    @test isstructmatch(iso, iso2, :exact, prefilter=false)
    @test !isstructmatch(cyclohexane, hexane, :exact, prefilter=false)
    @test !isstructmatch(cyclohexane, hexane, :nodeinduced, prefilter=false)
    @test isstructmatch(cyclohexane, hexane, :edgeinduced, prefilter=false)

    tms = smilestomol("C[Si](C)(C)C")
    tsm = smilestomol("[Si]C([Si])([Si])[Si]")
    @test !isstructmatch(tms, tsm, :exact, prefilter=false)

    sulfide = smilestomol("CSSC")
    disconn = smilestomol("CS.SC")
    @test !isstructmatch(sulfide, disconn, :exact, prefilter=false)
    @test isstructmatch(sulfide, disconn, :edgeinduced, prefilter=false)

    # Invalid (no edges)
    nacl = smilestomol("[Na+].[Cl-]")
    na = smilestomol("[Na]")
    @test !isstructmatch(nacl, na, :edgeinduced, prefilter=false)

    tetrahedrane = smilestomol("C12C3C1C23")
    fused = smilestomol("C1C2C1C2")
    spiro = smilestomol("C1CC12CC2")
    @test isstructmatch(tetrahedrane, fused, :edgeinduced, prefilter=false)
    @test !isstructmatch(tetrahedrane, spiro, :edgeinduced, prefilter=false)
end

@testset "connectedquery" begin
    hexane = smilestomol("CCCCCC")
    anyatom = parse(SMARTS, "*")
    @test isstructmatch(hexane, anyatom, :substruct)

    aniline = smilestomol("C1=CC=CC=C1N")
    dieamine = smilestomol("CCNCC")
    priamine = parse(SMARTS, "[NX3;H2]")
    @test isstructmatch(aniline, priamine, :substruct)
    @test !isstructmatch(dieamine, priamine, :substruct)

    diether = smilestomol("CCOCC")
    phenol = smilestomol("c1ccccc1O")
    glycerol = smilestomol("OCC(O)CO")
    alcohol = parse(SMARTS, "[#6][OD]")
    @test !isstructmatch(diether, alcohol, :substruct)
    @test isstructmatch(phenol, alcohol, :substruct)
    @test isstructmatch(glycerol, alcohol, :substruct)

    cyclopentane = smilestomol("C1CCCC1")
    pyrrole = smilestomol("N1C=CC=C1")
    aliphring = parse(SMARTS, "*@;!:*")
    @test isstructmatch(cyclopentane, aliphring, :substruct)
    @test !isstructmatch(pyrrole, aliphring, :substruct)

    triamine = smilestomol("CCN(CC)CC")
    acetamide = smilestomol("CC(=O)N")
    notamide = parse(SMARTS, raw"[NX3;!$(NC=O)]")
    @test isstructmatch(triamine, notamide, :substruct)
    @test !isstructmatch(acetamide, notamide, :substruct)

    po1 = smilestomol("COOC")
    npo1 = smilestomol("COCOCOC")
    peroxide = parse(SMARTS, "[OX2][OX2]")
    @test isstructmatch(po1, peroxide, :substruct)
    @test !isstructmatch(npo1, peroxide, :substruct)

    pyridine = smilestomol("n1ccccc1")
    pyrrole = smilestomol("[nH]1cccc1")
    sixmem = parse(SMARTS, "[*r6]1[*r6][*r6][*r6][*r6][*r6]1")
    @test isstructmatch(pyridine, sixmem, :substruct)
    @test !isstructmatch(pyrrole, sixmem, :substruct)
end

@testset "disconnectedquery" begin
    hetero1 = smilestomol("CCN(O)[O-]")
    hetero2 = smilestomol("CCN=N")
    hetero3 = smilestomol("BrCCOC#N")
    hetero4 = smilestomol("CCS(=O)(=O)O")
    disconn = parse(SMARTS, "[#7,#8].[!#6].N")
    @test isstructmatch(hetero1, disconn, :substruct)
    @test !isstructmatch(hetero2, disconn, :substruct)
    @test isstructmatch(hetero3, disconn, :substruct)
    @test !isstructmatch(hetero4, disconn, :substruct)
end

@testset "node matching" begin
    function nodematch(mol, query, idx=1, len=1)
        afunc = atommatch(mol, query)
        bfunc = bondmatch(mol, query)
        emaps = edgesubgraphmatches(mol, query, nodematcher=afunc, edgematcher=bfunc)
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

end # substructure
