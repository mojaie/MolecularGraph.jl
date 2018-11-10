
@testset "substructure" begin

@testset "filter" begin
    propane = smilestomol("CCC")
    butane = smilestomol("CCCC")
    @test !fast_identity_filter(propane, butane)
    @test fast_substr_filter(butane, propane)
    @test !fast_substr_filter(propane, butane)
    BuOH = smilestomol("CCCC")
    @test fast_identity_filter(butane, BuOH)
end

@testset "substruct" begin
    hexane = smilestomol("CCCCCC")
    iso = smilestomol("CC(C)CCC")
    iso2 = smilestomol("CCCC(C)C")
    @test !is_identical(hexane, iso)
    @test is_identical(iso, iso2)
    cyclohexane = smilestomol("C1CCCCC1")
    @test !is_identical(hexane, cyclohexane)
    @test is_substruct(hexane, cyclohexane)
    tms = smilestomol("C[Si](C)(C)C")
    tsm = smilestomol("[Si]C([Si])([Si])[Si]")
    @test !is_identical(tms, tsm)
    sulfide = smilestomol("CSSC")
    disconn = smilestomol("CS.SC")
    @test !is_identical(sulfide, disconn)
    @test is_substruct(disconn, sulfide)
end

@testset "deltaY" begin
    cyclopropane = smilestomol("C1CC1")
    isopropane = smilestomol("CC(C)C")
    @test !is_identical(cyclopropane, isopropane)
    tetrahedrane = smilestomol("C12C3C1C23")
    fused = smilestomol("C1C2C1C2")
    spiro = smilestomol("C1CC12CC2")
    @test is_substruct(fused, tetrahedrane)
    @test !is_substruct(spiro, tetrahedrane)
end

end # substructure
