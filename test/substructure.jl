
@testset "substructure" begin

@testset "filter" begin
    BuOH = smilestomol("CCCO")
    butane = smilestomol("CCCC")
    propane = smilestomol("CCC")
    @test fast_identity_filter(butane, BuOH)
    @test !fast_identity_filter(propane, butane)
    @test fast_substr_filter(butane, propane)
    @test !fast_substr_filter(propane, butane)
end

@testset "substruct" begin
    null = smilestomol("")
    @test !is_identical(null, null)
    @test !is_substruct(null, null)

    hexane = smilestomol("CCCCCC")
    iso = smilestomol("CC(C)CCC")
    iso2 = smilestomol("CCCC(C)C")
    cyclohexane = smilestomol("C1CCCCC1")
    @test !is_identical(hexane, iso)
    @test is_identical(iso, iso2)
    @test !is_identical(hexane, cyclohexane)
    @test is_substruct(hexane, cyclohexane)

    tms = smilestomol("C[Si](C)(C)C")
    tsm = smilestomol("[Si]C([Si])([Si])[Si]")
    @test !is_identical(tms, tsm)

    sulfide = smilestomol("CSSC")
    disconn = smilestomol("CS.SC")
    @test !is_identical(sulfide, disconn)
    @test is_substruct(disconn, sulfide)

    # Invalid (no edges)
    nacl = smilestomol("[Na+].[Cl-]")
    na = smilestomol("[Na]")
    @test !is_substruct(na, nacl)

    tetrahedrane = smilestomol("C12C3C1C23")
    fused = smilestomol("C1C2C1C2")
    spiro = smilestomol("C1CC12CC2")
    @test is_substruct(fused, tetrahedrane)
    @test !is_substruct(spiro, tetrahedrane)
end

@testset "connectedquery" begin
    anyatom = parse(ConnectedSMARTS, "*")
    hexane = smilestomol("CCCCCC")
    @test is_querymatch(hexane, anyatom)

    priamine = parse(ConnectedSMARTS, "[NX3;H2]")
    aniline = smilestomol("C1=CC=CC=C1N")
    dieamine = smilestomol("CCNCC")
    @test is_querymatch(aniline, priamine)
    @test !is_querymatch(dieamine, priamine)

    alcohol = parse(ConnectedSMARTS, "[#6][OD]")
    diether = smilestomol("CCOCC")
    phenol = smilestomol("c1ccccc1O")
    glycerol = smilestomol("OCC(O)CO")
    @test !is_querymatch(diether, alcohol)
    @test is_querymatch(phenol, alcohol)
    @test is_querymatch(glycerol, alcohol)

    aliphring = parse(ConnectedSMARTS, "*@;!:*")
    cyclopentane = smilestomol("C1CCCC1")
    pyrrole = smilestomol("N1C=CC=C1")
    @test is_querymatch(cyclopentane, aliphring)
    @test !is_querymatch(pyrrole, aliphring)

    notamide = parse(ConnectedSMARTS, raw"[NX3;!$(NC=O)]")
    triamine = smilestomol("CCN(CC)CC")
    acetamide = smilestomol("CC(=O)N")
    @test is_querymatch(triamine, notamide)
    @test !is_querymatch(acetamide, notamide)

    peroxide = parse(ConnectedSMARTS, "[OX2][OX2]")
    po1 = smilestomol("COOC")
    npo1 = smilestomol("COCOCOC")
    @test is_querymatch(po1, peroxide)
    @test !is_querymatch(npo1, peroxide)

    sixmem = parse(ConnectedSMARTS, "[*r6]1[*r6][*r6][*r6][*r6][*r6]1")
    pyridine = smilestomol("n1ccccc1")
    pyrrole = smilestomol("n1cccc1")
    @test is_querymatch(pyridine, sixmem)
    @test !is_querymatch(pyrrole, sixmem)
end

@testset "disconnectedquery" begin
    disconn = parse(SMARTS, "[#7,#8].[!#6].N")
    hetero1 = smilestomol("CCN(O)[O-]")
    hetero2 = smilestomol("CCN=N")
    hetero3 = smilestomol("BrCCOC#N")
    hetero4 = smilestomol("CCS(=O)(=O)O")
    @test is_querymatch(hetero1, disconn)
    @test !is_querymatch(hetero2, disconn)
    @test is_querymatch(hetero3, disconn)
    @test !is_querymatch(hetero4, disconn)
end

end # substructure
