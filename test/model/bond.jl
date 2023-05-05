
@testset "model.bond" begin

@testset "bond" begin

    sdfb = SDFBond(2, 0, true)
    @test sdfb[:order] == 2
    @test sdfb[:notation] == 0
    @test sdfb[:isordered]
    sdfb2 = SDFBond([2, 0, true])
    @test hash(sdfb) == hash(sdfb2)
    @test sdfb == sdfb2

    smib = SMILESBond(1, false, :up)
    @test smib[:order] == 1
    @test smib[:isaromatic] == 0
    @test smib[:direction] === :up
    smib2 = SMILESBond([1, false, "up"])
    @test hash(smib) == hash(smib2)
    @test smib == smib2
end

end # model.bond