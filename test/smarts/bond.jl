
using MolecularGraph: bondsymbol!, bond!


@testset "smarts.bond" begin

@testset "bondsymbol" begin
    state = SmartsParser("", false)
    implicit1 = bondsymbol!(state)
    @test implicit1 === nothing

    state = SmartsParser("-", false)
    explicit1 = bondsymbol!(state)
    @test explicit1 == QueryFormula(:and, Set([
        QueryFormula(:bondorder, 1),
        QueryFormula(:isaromaticbond, false)
    ]))

    state = SmartsParser("\\?", false)
    stereo4 = bondsymbol!(state)
    @test stereo4 == QueryFormula(:not, QueryFormula(:stereo, :up))
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
    @test anyb.query == QueryFormula(:any, true)

    state = SmartsParser("-!@", false)
    notring = bond!(state)
    @test notring.query == QueryFormula(:and, Set([
        QueryFormula(:bondorder, 1),
        QueryFormula(:isaromaticbond, false),
        QueryFormula(:isringbond, false)
    ]))
end

end # smiles.bond
