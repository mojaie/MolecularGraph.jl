
@testset "smarts.smarts" begin

@testset "chain" begin
    nullmol = parse(SMARTS, "")
    @test nodecount(nullmol) == 0
    @test edgecount(nullmol) == 0

    aliphatic = parse(SMARTS, "C")
    @test nodeattr(aliphatic, 1).query == QueryFormula(:and, Set([
        QueryFormula(:atomsymbol, :C),
        QueryFormula(:isaromatic, false)
    ]))
    @test edgecount(aliphatic) == 0

    carbonyl = parse(SMARTS, "[CX3]=[OX1]")
    @test nodeattr(carbonyl, 1).query == QueryFormula(:and, Set([
        QueryFormula(:atomsymbol, :C),
        QueryFormula(:isaromatic, false),
        QueryFormula(:connectivity, 3)
    ]))
    @test nodeattr(carbonyl, 2).query == QueryFormula(:and, Set([
        QueryFormula(:atomsymbol, :O),
        QueryFormula(:isaromatic, false),
        QueryFormula(:connectivity, 1)
    ]))

    ether = parse(SMARTS, "[#6][OD2][#6]")
    @test nodeattr(ether, 1).query == QueryFormula(:atomsymbol, :C)
    @test nodeattr(ether, 2).query == QueryFormula(:and, Set([
        QueryFormula(:atomsymbol, :O),
        QueryFormula(:isaromatic, false),
        QueryFormula(:nodedegree, 2)
    ]))

    notH = parse(SMARTS, "[!#1]")
    @test nodeattr(notH, 1).query == QueryFormula(:not, QueryFormula(:atomsymbol, :H))

    imine = parse(SMARTS, raw"[$([CX3]([#6])[#6]),$([CX3H][#6])]=[$([NX2][#6]),$([NX2H])]")
    @test nodeattr(imine, 1).query == QueryFormula(:or, Set([
        QueryFormula(:recursive, "[CX3]([#6])[#6]"),
        QueryFormula(:recursive, "[CX3H][#6]")
    ]))
    @test nodeattr(imine, 2).query == QueryFormula(:or, Set([
        QueryFormula(:recursive, "[NX2][#6]"),
        QueryFormula(:recursive, "[NX2H]")
    ]))
end

end # smarts.SMARTS
