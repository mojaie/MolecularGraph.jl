
@testset "mcs" begin

@testset "MCES" begin
    null = smilestomol("")
    @test mcssize(null, null) == 0

    BuOH = smilestomol("CCCO")
    butane = smilestomol("CCCC")
    @test mcssize(butane, BuOH) == 2

    hexane = smilestomol("CCCCCC")
    cyclohexane = smilestomol("C1CCCCC1")
    @test mcssize(hexane, cyclohexane) == 5

    tms = smilestomol("C[Si](C)(C)C")
    tsm = smilestomol("[Si]C([Si])([Si])[Si]")
    @test mcssize(tms, tsm) == 1

    sulfide = smilestomol("CSSC")
    disconn = smilestomol("CS.SC")
    @test mcssize(sulfide, disconn) == 2

    # Delta-Y
    cyclopropane = smilestomol("C1CC1")
    isopropane = smilestomol("CC(C)C")
    @test mcssize(cyclopropane, isopropane) == 2

    tetrahedrane = smilestomol("C12C3C1C23")
    fused = smilestomol("C1C2C1C2")
    spiro = smilestomol("C1CC12CC2")
    @test mcssize(tetrahedrane, fused) == 5
    @test mcssize(tetrahedrane, spiro) == 4
end

end # mcs
