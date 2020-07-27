
using MolecularGraph.Graph: lgnodematcher, lgedgematcher
using MolecularGraph:
    fastidentityfilter, fastsubstrfilter,
    atommatch, bondmatch

@testset "substructure" begin

@testset "filter" begin
    BuOH = smilestomol("CCCO")
    butane = smilestomol("CCCC")
    propane = smilestomol("CCC")
    @test fastidentityfilter(butane, BuOH)
    @test !fastidentityfilter(propane, butane)
    @test fastsubstrfilter(butane, propane)
    @test !fastsubstrfilter(propane, butane)
end


@testset "substruct" begin
    null = smilestomol("")
    @test !isstructmatch(null, null, prefilter=false)
    @test !issubstructmatch(null, null, prefilter=false)

    hexane = smilestomol("CCCCCC")
    iso = smilestomol("CC(C)CCC")
    iso2 = smilestomol("CCCC(C)C")
    cyclohexane = smilestomol("C1CCCCC1")
    @test !isstructmatch(hexane, iso, prefilter=false)
    @test isstructmatch(iso, iso2, prefilter=false)
    @test !isstructmatch(cyclohexane, hexane, prefilter=false)
    @test issubstructmatch(cyclohexane, hexane, prefilter=false)

    tms = smilestomol("C[Si](C)(C)C")
    tsm = smilestomol("[Si]C([Si])([Si])[Si]")
    @test !isstructmatch(tms, tsm, prefilter=false)

    sulfide = smilestomol("CSSC")
    disconn = smilestomol("CS.SC")
    @test !isstructmatch(sulfide, disconn, prefilter=false)
    @test issubstructmatch(sulfide, disconn, prefilter=false)

    # Invalid (no edges)
    nacl = smilestomol("[Na+].[Cl-]")
    na = smilestomol("[Na]")
    @test !issubstructmatch(nacl, na, prefilter=false)

    tetrahedrane = smilestomol("C12C3C1C23")
    fused = smilestomol("C1C2C1C2")
    spiro = smilestomol("C1CC12CC2")
    @test issubstructmatch(tetrahedrane, fused, prefilter=false)
    @test !issubstructmatch(tetrahedrane, spiro, prefilter=false)
end

@testset "connectedquery" begin
    anyatom = parse(SMARTS, "*")
    hexane = smilestomol("CCCCCC")
    @test isquerymatch(hexane, anyatom)

    priamine = parse(SMARTS, "[NX3;H2]")
    aniline = smilestomol("C1=CC=CC=C1N")
    dieamine = smilestomol("CCNCC")
    @test isquerymatch(aniline, priamine)
    @test !isquerymatch(dieamine, priamine)

    alcohol = parse(SMARTS, "[#6][OD]")
    diether = smilestomol("CCOCC")
    phenol = smilestomol("c1ccccc1O")
    glycerol = smilestomol("OCC(O)CO")
    @test !isquerymatch(diether, alcohol)
    @test isquerymatch(phenol, alcohol)
    @test isquerymatch(glycerol, alcohol)

    aliphring = parse(SMARTS, "*@;!:*")
    cyclopentane = smilestomol("C1CCCC1")
    pyrrole = smilestomol("N1C=CC=C1")
    @test isquerymatch(cyclopentane, aliphring)
    @test !isquerymatch(pyrrole, aliphring)

    notamide = parse(SMARTS, raw"[NX3;!$(NC=O)]")
    triamine = smilestomol("CCN(CC)CC")
    acetamide = smilestomol("CC(=O)N")
    @test isquerymatch(triamine, notamide)
    @test !isquerymatch(acetamide, notamide)

    peroxide = parse(SMARTS, "[OX2][OX2]")
    po1 = smilestomol("COOC")
    npo1 = smilestomol("COCOCOC")
    @test isquerymatch(po1, peroxide)
    @test !isquerymatch(npo1, peroxide)

    sixmem = parse(SMARTS, "[*r6]1[*r6][*r6][*r6][*r6][*r6]1")
    pyridine = smilestomol("n1ccccc1")
    pyrrole = smilestomol("[nH]1cccc1")
    @test isquerymatch(pyridine, sixmem)
    @test !isquerymatch(pyrrole, sixmem)
end

@testset "disconnectedquery" begin
    disconn = parse(SMARTS, "[#7,#8].[!#6].N")
    hetero1 = smilestomol("CCN(O)[O-]")
    hetero2 = smilestomol("CCN=N")
    hetero3 = smilestomol("BrCCOC#N")
    hetero4 = smilestomol("CCS(=O)(=O)O")
    @test isquerymatch(hetero1, disconn)
    @test !isquerymatch(hetero2, disconn)
    @test isquerymatch(hetero3, disconn)
    @test !isquerymatch(hetero4, disconn)
end

end # substructure
