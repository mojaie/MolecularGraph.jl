

@testset "structurematch_mcs" begin

@testset "MCS" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    null = smilestomol("")
    @test length(disconnected_mcis(null, null)[1]) == 0
    @test length(disconnected_mces(null, null)[1]) == 0
    @test length(connected_mcis(null, null)[1]) == 0
    @test length(connected_mces(null, null)[1]) == 0
    @test length(tdmcis(null, null)[1]) == 0
    @test length(tdmcis(null, null, tolerance=1)[1]) == 0
    @test length(tdmces(null, null)[1]) == 0
    @test length(tdmces(null, null, tolerance=1)[1]) == 0

    BuOH = smilestomol("CCCO")
    butane = smilestomol("CCCC")
    @test length(disconnected_mcis(butane, BuOH)[1]) == 3
    @test length(disconnected_mces(butane, BuOH)[1]) == 2

    hexane = smilestomol("CCCCCC")
    cyclohexane = smilestomol("C1CCCCC1")
    @test length(tdmcis(hexane, cyclohexane)[1]) == 4
    @test length(tdmcis(hexane, cyclohexane, tolerance=1)[1]) == 4
    @test length(tdmces(hexane, cyclohexane)[1]) == 4
    @test length(tdmces(hexane, cyclohexane, tolerance=1)[1]) == 4

    tms = smilestomol("C[Si](C)(C)C")
    tsm = smilestomol("[Si]C([Si])([Si])[Si]")
    @test length(disconnected_mcis(tms, tsm)[1]) == 2
    @test length(disconnected_mces(tms, tsm)[1]) == 1
    @test length(tdmcis(tms, tsm)[1]) == 2
    @test length(tdmces(tms, tsm)[1]) == 1

    sulfide = smilestomol("CSSC")
    disconn = smilestomol("CS.SC")
    @test length(connected_mcis(sulfide, disconn)[1]) == 2
    @test length(disconnected_mcis(sulfide, disconn)[1]) == 3
    @test length(connected_mces(sulfide, disconn)[1]) == 1
    @test length(disconnected_mces(sulfide, disconn)[1]) == 2
    @test length(tdmcis(sulfide, disconn, tolerance=1)[1]) == 2
    @test length(tdmces(sulfide, disconn, tolerance=1)[1]) == 1

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
    @test length(tdmcis(cid6437877, cid5481173, tolerance=1)[1]) == 28
    @test length(tdmces(cid6437877, cid5481173, tolerance=1)[1]) == 29
    @test length(tdmcis(cid6437877, cid5481173, diameter=8)[1]) == 18
    @test length(tdmces(cid6437877, cid5481173, diameter=8)[1]) == 20
    # global_logger(default_logger)
end

end # structurematch_mcs
