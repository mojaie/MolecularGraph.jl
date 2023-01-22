
using MolecularGraph: bondsymbol!, bond!


@testset "smarts.bond" begin

@testset "bondsymbol" begin
    state = SMARTSParser("")
    implicit1 = bondsymbol!(state)
    @test implicit1 === nothing

    state = SMARTSParser("-")
    explicit1 = QueryTruthTable(bondsymbol!(state)...)
    @test explicit1 == QueryTruthTable(v -> v[1] & ~v[2], [(:order, 1), (:isaromatic,)])

    state = SMARTSParser(raw"\?")
    stereo4 = QueryTruthTable(bondsymbol!(state)...)
    @test stereo4 == QueryTruthTable(v -> ~v[1], [(:stereo, :up)])
    @test state.pos == 3
end

@testset "bond" begin
    state = SMILESParser("#")
    triple = bond!(state)
    @test triple[:order] == 3

    state = SMILESParser(":")
    arom = bond!(state)
    @test arom[:isaromatic]
end

@testset "smartsbond" begin
    state = SMARTSParser("~")
    anyb = bond!(state)
    @test anyb == QueryTruthTable(any_query(true)...)

    state = SMARTSParser("-!@")
    notring = bond!(state)
    @test notring == QueryTruthTable(
        v -> v[1] & ~v[2] & ~v[3],
        [(:order, 1), (:isaromatic,), (:isring,)])
end

end # smiles.bond
