

@testset "structurematch_mcs" begin

@testset "MCS" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    null = smilestomol("")
    @test length(disconnected_mcis(null, null)[1]) == 0

    BuOH = smilestomol("CCCO")
    butane = smilestomol("CCCC")
    @test length(disconnected_mcis(butane, BuOH)[1]) == 3
    @test length(disconnected_mces(butane, BuOH)[1]) == 2

    hexane = smilestomol("CCCCCC")
    cyclohexane = smilestomol("C1CCCCC1")
    @test length(disconnected_mcis(hexane, cyclohexane)[1]) == 5
    @test length(disconnected_mces(hexane, cyclohexane)[1]) == 5

    tms = smilestomol("C[Si](C)(C)C")
    tsm = smilestomol("[Si]C([Si])([Si])[Si]")
    @test length(disconnected_mcis(tms, tsm)[1]) == 2
    @test length(disconnected_mces(tms, tsm)[1]) == 1

    sulfide = smilestomol("CSSC")
    disconn = smilestomol("CS.SC")
    @test length(connected_mcis(sulfide, disconn)[1]) == 2
    @test length(connected_mces(sulfide, disconn)[1]) == 1
    @test length(disconnected_mcis(sulfide, disconn)[1]) == 3
    @test length(disconnected_mces(sulfide, disconn)[1]) == 2
    @test length(tdmcis(sulfide, disconn)[1]) == 2
    @test length(tdmces(sulfide, disconn)[1]) == 1

    # Delta-Y
    cyclopropane = smilestomol("C1CC1")
    isopropane = smilestomol("CC(C)C")
    @test length(disconnected_mcis(cyclopropane, isopropane)[1]) == 2
    @test length(disconnected_mces(cyclopropane, isopropane)[1]) == 2

    tetrahedrane = smilestomol("C12C3C1C23")
    fused = smilestomol("C1C2C1C2")
    spiro = smilestomol("C1CC12CC2")
    @test length(disconnected_mcis(tetrahedrane, fused)[1]) == 3
    @test length(disconnected_mcis(tetrahedrane, spiro)[1]) == 3
    @test length(disconnected_mces(tetrahedrane, fused)[1]) == 5
    @test length(disconnected_mces(tetrahedrane, spiro)[1]) == 4

    cid6437877 = smilestomol("CC1=C(SC=N1)C=CC2=C(N3C(C(C3=O)NC(=O)C(=NOC)C4=CSC(=N4)N)SC2)C(=O)OCOC(=O)C(C)(C)C")  # Cefditoren pivoxil
    cid5481173 = smilestomol("CC(C)(C(=O)O)ON=C(C1=CSC(=N1)N)C(=O)NC2C3N(C2=O)C(=C(CS3)C[N+]4=CC=CC=C4)C(=O)[O-]")  # Ceftazidime
    @test length(tdmcis(cid6437877, cid5481173, tolerance=1)[1]) == 27
    @test length(tdmces(cid6437877, cid5481173, tolerance=1)[1]) == 28
    @test length(tdmcis(cid6437877, cid5481173, diameter=8)[1]) == 18
    @test length(tdmces(cid6437877, cid5481173, diameter=8)[1]) == 20
    # global_logger(default_logger)
end

@testset "MCS constraint array" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    null = smilestomol("")
    nullc = mcis_constraints(null)
    @test length(maximum_common_subgraph(nullc, nullc)[1]) == 0
    nullc = mces_constraints(null)
    @test length(maximum_common_subgraph(nullc, nullc)[1]) == 0
    nullc = tdmcis_constraints(null)
    @test length(maximum_common_subgraph(nullc, nullc)[1]) == 0
    @test length(maximum_common_subgraph(nullc, nullc, tolerance=1)[1]) == 0
    nullc = tdmces_constraints(null)
    @test length(maximum_common_subgraph(nullc, nullc)[1]) == 0
    @test length(maximum_common_subgraph(nullc, nullc, tolerance=1)[1]) == 0
    nullc = tdmcis_constraints(null, :any)
    @test length(maximum_common_subgraph(nullc, nullc)[1]) == 0
    nullc = tdmces_constraints(null, :any)
    @test length(maximum_common_subgraph(nullc, nullc)[1]) == 0

    BuOH = smilestomol("CCCO")
    butane = smilestomol("CCCC")
    BuOHc = mcis_constraints(BuOH)
    butanec = mcis_constraints(butane)
    @test length(maximum_common_subgraph(butanec, BuOHc)[1]) == 3
    BuOHc = mces_constraints(BuOH)
    butanec = mces_constraints(butane)
    @test length(maximum_common_subgraph(butanec, BuOHc)[1]) == 2

    hexane = smilestomol("CCCCCC")
    cyclohexane = smilestomol("C1CCCCC1")
    hexanec = tdmcis_constraints(hexane)
    cyclohexanec = tdmcis_constraints(cyclohexane)
    @test length(maximum_common_subgraph(hexanec, cyclohexanec)[1]) == 4
    @test length(maximum_common_subgraph(hexanec, cyclohexanec, tolerance=1)[1]) == 4
    hexanec = tdmces_constraints(hexane)
    cyclohexanec = tdmces_constraints(cyclohexane)
    @test length(maximum_common_subgraph(hexanec, cyclohexanec)[1]) == 4
    @test length(maximum_common_subgraph(hexanec, cyclohexanec, tolerance=1)[1]) == 4
    hexanec = tdmcis_constraints(hexane, :any)
    cyclohexanec = tdmcis_constraints(cyclohexane, :any)
    @test length(maximum_common_subgraph(hexanec, cyclohexanec)[1]) == 5
    hexanec = tdmces_constraints(hexane, :any)
    cyclohexanec = tdmces_constraints(cyclohexane, :any)
    @test length(maximum_common_subgraph(hexanec, cyclohexanec)[1]) == 5

    tms = smilestomol("C[Si](C)(C)C")
    tsm = smilestomol("[Si]C([Si])([Si])[Si]")
    tmsc = mcis_constraints(tms)
    tsmc = mcis_constraints(tsm)
    @test length(maximum_common_subgraph(tmsc, tsmc)[1]) == 2
    tmsc = mces_constraints(tms)
    tsmc = mces_constraints(tsm)
    @test length(maximum_common_subgraph(tmsc, tsmc)[1]) == 1
    tmsc = tdmcis_constraints(tms)
    tsmc = tdmcis_constraints(tsm)
    @test length(maximum_common_subgraph(tmsc, tsmc)[1]) == 2
    tmsc = tdmces_constraints(tms)
    tsmc = tdmces_constraints(tsm)
    @test length(maximum_common_subgraph(tmsc, tsmc)[1]) == 1

    sulfide = smilestomol("CSSC")
    disconn = smilestomol("CS.SC")
    sulfidec = mcis_constraints(sulfide)
    disconnc = mcis_constraints(disconn)
    # @test length(maximum_common_subgraph(sulfidec, disconnc, connected=true)[1]) == 2
    @test length(maximum_common_subgraph(sulfidec, disconnc)[1]) == 3
    sulfidec = mces_constraints(sulfide)
    disconnc = mces_constraints(disconn)
    # @test length(maximum_common_edge_subgraph(sulfidec, disconnc, connected=true)[1]) == 1
    @test length(maximum_common_subgraph(sulfidec, disconnc)[1]) == 2
    sulfidec = tdmcis_constraints(sulfide)
    disconnc = tdmcis_constraints(disconn)
    @test length(maximum_common_subgraph(sulfidec, disconnc, tolerance=1)[1]) == 2
    sulfidec = tdmces_constraints(sulfide)
    disconnc = tdmces_constraints(disconn)
    @test length(maximum_common_subgraph(sulfidec, disconnc, tolerance=1)[1]) == 1

    # Delta-Y
    cyclopropane = smilestomol("C1CC1")
    isopropane = smilestomol("CC(C)C")
    cyclopropanec = mcis_constraints(cyclopropane)
    isopropanec = mcis_constraints(isopropane)
    @test length(maximum_common_subgraph(cyclopropanec, isopropanec)[1]) == 2
    cyclopropanec = mces_constraints(cyclopropane)
    isopropanec = mces_constraints(isopropane)
    @test length(maximum_common_subgraph(cyclopropanec, isopropanec)[1]) == 2

    tetrahedrane = smilestomol("C12C3C1C23")
    fused = smilestomol("C1C2C1C2")
    spiro = smilestomol("C1CC12CC2")
    tetrahedranec = mcis_constraints(tetrahedrane)
    fusedc = mcis_constraints(fused)
    spiroc = mcis_constraints(spiro)
    @test length(maximum_common_subgraph(tetrahedranec, fusedc)[1]) == 3
    @test length(maximum_common_subgraph(tetrahedranec, spiroc)[1]) == 3
    tetrahedranec = mces_constraints(tetrahedrane)
    fusedc = mces_constraints(fused)
    spiroc = mces_constraints(spiro)
    @test length(maximum_common_subgraph(tetrahedranec, fusedc)[1]) == 5
    @test length(maximum_common_subgraph(tetrahedranec, spiroc)[1]) == 4

    cid6437877 = smilestomol("CC1=C(SC=N1)C=CC2=C(N3C(C(C3=O)NC(=O)C(=NOC)C4=CSC(=N4)N)SC2)C(=O)OCOC(=O)C(C)(C)C")  # Cefditoren pivoxil
    cid5481173 = smilestomol("CC(C)(C(=O)O)ON=C(C1=CSC(=N1)N)C(=O)NC2C3N(C2=O)C(=C(CS3)C[N+]4=CC=CC=C4)C(=O)[O-]")  # Ceftazidime
    cid6437877t = tdmcis_constraints(cid6437877)
    cid5481173t = tdmcis_constraints(cid5481173)
    @test length(maximum_common_subgraph(cid6437877t, cid5481173t, tolerance=1)[1]) == 27
    cid6437877t = tdmces_constraints(cid6437877)
    cid5481173t = tdmces_constraints(cid5481173)
    @test length(maximum_common_subgraph(cid6437877t, cid5481173t, tolerance=1)[1]) == 28
    cid6437877t = tdmcis_constraints(cid6437877, diameter=8)
    cid5481173t = tdmcis_constraints(cid5481173, diameter=8)
    @test length(maximum_common_subgraph(cid6437877t, cid5481173t)[1]) == 18
    cid6437877t = tdmces_constraints(cid6437877, diameter=8)
    cid5481173t = tdmces_constraints(cid5481173, diameter=8)
    @test length(maximum_common_subgraph(cid6437877t, cid5481173t)[1]) == 20
    # global_logger(default_logger)
end

end # structurematch
