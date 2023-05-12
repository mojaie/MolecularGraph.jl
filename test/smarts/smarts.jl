
@testset "smarts.smarts" begin

@testset "smarts" begin
    SMARTSTT = MolGraph{Int,QueryTruthTable,QueryTruthTable}
    nullmol = smartstomol("")
    @test nv(nullmol) == 0
    @test ne(nullmol) == 0

    aliphatic = smartstomol(SMARTSTT, "C")
    @test props(aliphatic, 1) == QueryTruthTable(
        v -> v[2] & ~v[1], [(:isaromatic,), (:symbol, :C)])
    @test ne(aliphatic) == 0

    carbonyl = smartstomol(SMARTSTT, "[CX3]=[OX1]")
    @test props(carbonyl, 1) == QueryTruthTable(
        v -> v[3] & ~v[2] & v[1], [(:connectivity, 3), (:isaromatic,), (:symbol, :C)])
    @test props(carbonyl, 2) == QueryTruthTable(
        v -> v[3] & ~v[2] & v[1], [(:connectivity, 1), (:isaromatic,), (:symbol, :O)])

    ether = smartstomol(SMARTSTT, "[#6][OD2][#6]")
    @test props(ether, 1) == QueryTruthTable(v -> v[1], [(:symbol, :C)])
    @test props(ether, 2) == QueryTruthTable(
        v -> v[3] & ~v[2] & v[1], [(:degree, 2), (:isaromatic,), (:symbol, :O)])

    notH = smartstomol(SMARTSTT, "[!#1]")
    @test props(notH, 1) == QueryTruthTable(v -> ~v[1], [(:symbol, :H)])

    imine = smartstomol(SMARTSTT, raw"[$([CX3]([#6])[#6]),$([CX3H][#6])]=[$([NX2][#6]),$([NX2H])]")
    @test props(imine, 1) == QueryTruthTable(
        v -> v[1] | v[2], [(:recursive, "[CX3H][#6]"), (:recursive, "[CX3]([#6])[#6]")])

    @test props(imine, 2) == QueryTruthTable(
        v -> v[1] | v[2], [(:recursive, "[NX2H]"), (:recursive, "[NX2][#6]")])
end

@testset "serialization" begin
    mol = smartstomol(raw"[$([CX3]([#6])[#6]),$([CX3H][#6])]=[$([NX2][#6]),$([NX2H])]")
    mol2 = SMARTSMolGraph(to_json(mol))
    @test mol == mol2
    @test mol !== mol2

    mol = smartstomol("[O,S]=P([O,S])([O,S])[O,S]")
    mol2 = SMARTSMolGraph(to_json(mol))
    @test mol == mol2
    @test mol !== mol2
end


end
