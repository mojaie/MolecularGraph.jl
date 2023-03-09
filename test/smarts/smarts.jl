
@testset "smarts.smarts" begin

@testset "chain" begin
    nullmol = smartstomol("")
    @test nv(nullmol) == 0
    @test ne(nullmol) == 0

    aliphatic = smartstomol("C")
    @test props(aliphatic, 1) == QueryTruthTable(
        v -> v[1] & ~v[2], [(:symbol, :C), (:isaromatic,)])
    @test ne(aliphatic) == 0

    carbonyl = smartstomol("[CX3]=[OX1]")
    @test props(carbonyl, 1) == QueryTruthTable(
        v -> v[1] & ~v[2] & v[3], [(:symbol, :C), (:isaromatic,), (:connectivity, 3)])
    @test props(carbonyl, 2) == QueryTruthTable(
        v -> v[1] & ~v[2] & v[3], [(:symbol, :O), (:isaromatic,), (:connectivity, 1)])

    ether = smartstomol("[#6][OD2][#6]")
    @test props(ether, 1) == QueryTruthTable(v -> v[1], [(:symbol, :C)])
    @test props(ether, 2) == QueryTruthTable(
        v -> v[1] & ~v[2] & v[3], [(:symbol, :O), (:isaromatic,), (:degree, 2)])

    notH = smartstomol("[!#1]")
    @test props(notH, 1) == QueryTruthTable(v -> ~v[1], [(:symbol, :H)])

    imine = smartstomol(raw"[$([CX3]([#6])[#6]),$([CX3H][#6])]=[$([NX2][#6]),$([NX2H])]")
    @test props(imine, 1) == QueryTruthTable(
        v -> v[1] | v[2], [(:recursive, "[CX3]([#6])[#6]"), (:recursive, "[CX3H][#6]")])

    @test props(imine, 2) == QueryTruthTable(
        v -> v[1] | v[2], [(:recursive, "[NX2][#6]"), (:recursive, "[NX2H]")])
end

end # smarts.SMARTS
