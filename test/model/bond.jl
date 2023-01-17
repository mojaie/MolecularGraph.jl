
@testset "model.bond" begin

@testset "bond" begin

    sdfb = SDFBond(2, 0, true)
    @test sdfb[:order] == 2
    @test sdfb[:notation] == 0
    @test sdfb[:is_ordered]

    smib = SMILESBond(1, false, :up)
    @test smib[:order] == 1
    @test smib[:is_aromatic] == 0
    @test smib[:direction] === :up
end

end # model.bond