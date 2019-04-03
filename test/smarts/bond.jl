
using MolecularGraph: bondsymbol!, bond!


@testset "smarts.bond" begin

@testset "bondsymbol" begin
    state = SmartsParser("", false)
    implicit1 = bondsymbol!(state)
    @test implicit1 === nothing

    state = SmartsParser("-", false)
    explicit1 = bondsymbol!(state)
    @test explicit1 == (:bondorder => 1)

    state = SmartsParser("\\?", false)
    stereo4 = bondsymbol!(state)
    @test stereo4 == (:stereo => 4)
    @test state.pos == 3
end

@testset "bond" begin
    state = SmilesParser("#", true)
    triple = bond!(state)
    @test triple.order == 3

    state = SmilesParser(":", true)
    arom = bond!(state)
    @test arom.isaromatic == true
end

@testset "smartsbond" begin
    state = SmartsParser("~", false)
    anyb = bond!(state)
    @test anyb.query == (:any => true)

    state = SmartsParser("-!@", false)
    notring = bond!(state)
    @test notring.query == (
        :and => (:bondorder => 1, :not => (:isringbond => true))
    )
end

end # smiles.bond
