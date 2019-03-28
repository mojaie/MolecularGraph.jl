
@testset "smarts.bond" begin

@testset "bondsymbol" begin
    state = SmartsParser{SMARTS}("", false)
    implicit1 = bondsymbol!(state)
    @test implicit1 === nothing

    state = SmartsParser{SMARTS}("-", false)
    explicit1 = bondsymbol!(state)
    @test explicit1 == (:bondorder => 1)

    state = SmartsParser{SMARTS}("\\?", false)
    stereo4 = bondsymbol!(state)
    @test stereo4 == (:stereo => 4)
    @test state.pos == 3
end

@testset "bond" begin
    state = SmartsParser{SMILES}("#", true)
    triple = bond!(state)
    @test triple.order == 3

    state = SmartsParser{SMILES}(":", true)
    arom = bond!(state)
    @test arom.isaromatic == true
end

@testset "smartsbond" begin
    state = SmartsParser{SMARTS}("~", false)
    anyb = bond!(state)
    @test anyb.query == (:any => true)

    state = SmartsParser{SMARTS}("-!@", false)
    notring = bond!(state)
    @test notring.query == (
        :and => (:bondorder => 1, :not => (:bond_isringmem => true))
    )
end

end # smiles.bond
