
@testset "mcs" begin

@testset "MCES" begin
    null = smilestomol("")
    @test mcesmolsize(null, null) == 0

    BuOH = smilestomol("CCCO")
    butane = smilestomol("CCCC")
    @test mcesmolsize(butane, BuOH) == 2

    hexane = smilestomol("CCCCCC")
    cyclohexane = smilestomol("C1CCCCC1")
    @test mcesmolsize(hexane, cyclohexane) == 5

    tms = smilestomol("C[Si](C)(C)C")
    tsm = smilestomol("[Si]C([Si])([Si])[Si]")
    @test mcesmolsize(tms, tsm) == 1

    sulfide = smilestomol("CSSC")
    disconn = smilestomol("CS.SC")
    @test mcesmolsize(sulfide, disconn) == 2

    # Delta-Y
    cyclopropane = smilestomol("C1CC1")
    isopropane = smilestomol("CC(C)C")
    @test mcesmolsize(cyclopropane, isopropane) == 2

    tetrahedrane = smilestomol("C12C3C1C23")
    fused = smilestomol("C1C2C1C2")
    spiro = smilestomol("C1CC12CC2")
    @test mcesmolsize(tetrahedrane, fused) == 5
    @test mcesmolsize(tetrahedrane, spiro) == 4
end

end # mcs
