
@testset "smarts.smarts" begin

@testset "smarts" begin
    nullmol = smartstomol("")
    @test nv(nullmol) == 0
    @test ne(nullmol) == 0

    aliphatic = smartstomol("C")
    aliphatic[1] == QueryAtom(
        [(1, 2), (1, 3), (3, 4)],
        [qand(), qeq(:symbol, "C"), qnot(), qtrue(:isaromatic)])
    @test ne(aliphatic) == 0

    carbonyl = smartstomol("[CX3]=[OX1]")
    @test carbonyl[1] == QueryAtom(
        [(1, 2), (1, 3), (3, 4), (3, 5), (5, 6)],
        [qand(), qeq(:connectivity, "3"), qand(), qeq(:symbol, "C"),
        qnot(), qtrue(:isaromatic)])
    @test carbonyl[2] == QueryAtom(
        [(1, 2), (1, 3), (3, 4), (3, 5), (5, 6)],
        [qand(), qeq(:connectivity, "1"), qand(), qeq(:symbol, "O"),
        qnot(), qtrue(:isaromatic)])

    ether = smartstomol("[#6][OD2][#6]")
    @test ether[1] == QueryAtom(Tuple{Int,Int}[], [qeq(:symbol, "C")])
    @test ether[2] == QueryAtom(
        [(1, 2), (1, 3), (3, 4), (3, 5), (5, 6)],
        [qand(), qeq(:degree, "2"), qand(), qeq(:symbol, "O"),
        qnot(), qtrue(:isaromatic)])

    notH = smartstomol("[!#1]")  # [!#1] -> [*] in initialization process
    @test notH[1] == QueryAtom(Tuple{Int,Int}[], [qanytrue()])

    imine = smartstomol(
        raw"[$([CX3]([#6])[#6]),$([CX3H][#6])]=[$([NX2][#6]),$([NX2H])]")
    @test imine[1] == QueryAtom(
        [(1, 2), (1, 3)],
        [qor(), qeq(:recursive, "[CX3H][#6]"), qeq(:recursive, "[CX3]([#6])[#6]")])

    @test imine[2] == QueryAtom(
        [(1, 2), (1, 3)],
        [qor(), qeq(:recursive, "[NX2H]"), qeq(:recursive, "[NX2][#6]")])
end

end
