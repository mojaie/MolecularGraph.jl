
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

    # TODO: substr match does not care charges
    nacl = smilestomol("[Na+].[Cl-]")
    na = smilestomol("[Na]")
    @test !is_identical(nacl, na)
    @test is_substruct(na, nacl)

    tetrahedrane = smilestomol("C12C3C1C23")
    fused = smilestomol("C1C2C1C2")
    spiro = smilestomol("C1CC12CC2")
    @test is_substruct(fused, tetrahedrane)
    @test !is_substruct(spiro, tetrahedrane)
end

@testset "query" begin
    anyatom = parse(SMARTS, "*")
    hexane = smilestomol("CCCCCC")
    @test match_molquery(hexane, anyatom)

    priamine = parse(SMARTS, "[NX3;H2]")
    aniline = smilestomol("C1=CC=CC=C1N")
    dieamine = smilestomol("CCNCC")
    @test match_molquery(aniline, priamine)
    @test !match_molquery(dieamine, priamine)

    alcohol = parse(SMARTS, "[#6][OD]")
    diether = smilestomol("CCOCC")
    phenol = smilestomol("c1ccccc1O")
    glycerol = smilestomol("OCC(O)CO")
    @test !match_molquery(diether, alcohol)
    @test match_molquery(phenol, alcohol)
    @test match_molquery(glycerol, alcohol)

    aliphring = parse(SMARTS, "*@;!:*")
    cyclopentane = smilestomol("C1CCCC1")
    pyrrole = smilestomol("N1C=CC=C1")
    @test match_molquery(cyclopentane, aliphring)
    @test !match_molquery(pyrrole, aliphring)

    notamide = parse(SMARTS, raw"[NX3;!$(NC=O)]")
    triamine = smilestomol("CCN(CC)CC")
    acetamide = smilestomol("CC(=O)N")
    @test match_molquery(triamine, notamide)
    @test !match_molquery(acetamide, notamide)

    disconn = parse(SMARTS, "[#7,#8].[!#6].N")
    hetero1 = smilestomol("CCN(O)[O-]")
    hetero2 = smilestomol("CCN=N")
    hetero3 = smilestomol("BrCCOC#N")
    hetero4 = smilestomol("CCS(=O)(=O)O")
    @test match_molquery(hetero1, disconn)
    @test !match_molquery(hetero2, disconn)
    @test match_molquery(hetero3, disconn)
    @test !match_molquery(hetero4, disconn)

    peroxide = parse(SMARTS, "[OX2][OX2]")
    po1 = smilestomol("COOC")
    npo1 = smilestomol("COCOCOC")
    @test match_molquery(po1, peroxide)
    @test !match_molquery(npo1, peroxide)
end

end # substructure
