
@testset "smarts.bond" begin

@testset "bondsymbol" begin
    state = SmartsParserState("")
    implicit1 = bondsymbol!(state)
    @test implicit1 === nothing

    state = SmartsParserState("-")
    explicit1 = bondsymbol!(state)
    @test explicit1 == (:BondOrder => 1)

    state = SmartsParserState("\\?")
    stereo4 = bondsymbol!(state)
    @test stereo4 == (:stereo => 4)
    @test state.pos == 3
end

@testset "bond" begin
    state = SmilesParserState("#")
    triple = bond!(state)
    @test triple.order == 3

    state = SmilesParserState(":")
    arom = bond!(state)
    @test arom.isaromatic == true
end

@testset "bondquery" begin
    state = SmartsParserState("~")
    anyb = bondquery!(state)
    @test anyb.query == (:any => true)

    state = SmartsParserState("-!@")
    notring = bondquery!(state)
    @test notring.query == (
        :and => (:BondOrder => 1, :not => (:RingBond => true))
    )
end

end # smiles.bond
