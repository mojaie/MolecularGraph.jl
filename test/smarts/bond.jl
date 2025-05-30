
@testset "smarts.bond" begin

@testset "bondsymbol" begin
    state = SMARTSParser{SMARTSMolGraph}("")
    implicit1 = bondsymbol!(state)
    @test implicit1 === nothing

    state = SMARTSParser{SMARTSMolGraph}("-")
    explicit1 = QueryTruthTable(bondsymbol!(state))
    @test explicit1 == QueryTruthTable(v -> v[2] & ~v[1], [(:isaromatic,), (:order, 1)])

    state = SMARTSParser{SMARTSMolGraph}(raw"\?")
    stereo4 = QueryTruthTable(bondsymbol!(state))
    @test stereo4 == QueryTruthTable(v -> ~v[1], [(:stereo, :up)])
    @test state.pos == 3
end

@testset "bond" begin
    state = SMILESParser{SMILESMolGraph}("#")
    triple = bond!(state)
    @test triple[:order] == 3

    state = SMILESParser{SMILESMolGraph}(":")
    arom = bond!(state)
    @test arom[:isaromatic]

    state = SMILESParser{SMILESMolGraph}("\\")
    down = bond!(state)
    @test down[:direction] === :down
end

@testset "smartsbond" begin
    SMARTSTT = MolGraph{Int,QueryTruthTable,QueryTruthTable}
    state = SMARTSParser{SMARTSTT}("~")
    anyb = bond!(state)
    @test anyb == QueryTruthTable(v -> true, [])

    state = SMARTSParser{SMARTSTT}("-!@")
    notring = bond!(state)
    @test notring == QueryTruthTable(
        v -> v[3] & ~v[1] & ~v[2], [(:is_in_ring,), (:isaromatic,), (:order, 1)])
end

end # smiles.bond
