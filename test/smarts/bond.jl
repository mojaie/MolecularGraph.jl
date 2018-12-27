
@testset "smarts.bond" begin

@testset "bondsymbol" begin
    state = ConnectedSmarts("")
    implicit1 = bondsymbol!(state)
    @test implicit1 === nothing

    state = ConnectedSmarts("-")
    explicit1 = bondsymbol!(state)
    @test explicit1 == (:BondOrder => 1)

    state = ConnectedSmarts("\\?")
    stereo4 = bondsymbol!(state)
    @test stereo4 == (:stereo => 4)
    @test state.pos == 3
end

@testset "bond" begin
    state = SmilesParser("#")
    triple = bond!(state)
    @test triple.order == 3

    state = SmilesParser(":")
    arom = bond!(state)
    @test arom.isaromatic == true
end

@testset "smartsbond" begin
    state = ConnectedSmarts("~")
    anyb = bond!(state)
    @test anyb.query == (:any => true)

    state = ConnectedSmarts("-!@")
    notring = bond!(state)
    @test notring.query == (
        :and => (:BondOrder => 1, :not => (:RingBond => true))
    )
end

end # smiles.bond
