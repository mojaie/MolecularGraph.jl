
@testset "model.bond" begin

@testset "bond" begin

    sdfb = SDFBond(; order=2, isordered=true)
    @test sdfb[:order] == 2
    @test sdfb[:notation] == 0
    @test sdfb[:isordered]
    sdfb2 = SDFBond(;order=2, isordered=true)
    @test sdfb == sdfb2
    @test hash(sdfb) == hash(sdfb2)

    smib = SMILESBond(; isaromatic=false, direction=:up)
    @test smib[:order] == 1
    @test smib[:isaromatic] == 0
    @test smib[:direction] === :up
    smib2 = SMILESBond(; isaromatic=false, direction=:up)
    @test smib == smib2
    @test hash(smib) == hash(smib2)
end

end # model.bond