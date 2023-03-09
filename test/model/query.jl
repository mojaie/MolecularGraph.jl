
@testset "model.query" begin

@testset "generation" begin
    gen1 = QueryTruthTable(
        and_query([
            statement(:symbol, :C),
            not_query(statement(:isaromatic)),
            or_query([
                statement(:connectivity, 3),
                statement(:connectivity, 4)
            ])
        ])...
    )
    func1 = QueryTruthTable(
        v -> v[1] & ~v[2] & (v[3] | v[4]),
        [(:symbol, :C), (:isaromatic,), (:connectivity, 3), (:connectivity, 4)]
    )
    @test issubset(gen1, func1)
    @test issubset(func1, gen1)
end

@testset "match" begin
    symC = QueryTruthTable(v -> v[1], [(:symbol, :C)])
    symN = QueryTruthTable(v -> v[1], [(:symbol, :N)])
    hc1 = QueryTruthTable(v -> v[1], [(:total_hydrogens, 1)])
    ch1 = QueryTruthTable(v -> v[1], [(:charge, 1)])
    @test issubset(symC, symC)
    @test !issubset(symC, symN)
    @test !issubset(symN, symC)
    @test issubset(hc1, hc1)
    @test !issubset(hc1, ch1)
    @test !issubset(ch1, hc1)

    notN = QueryTruthTable(v -> ~v[1], [(:symbol, :N)])
    @test !issubset(symN, notN)
    @test !issubset(notN, symN)

    and1 = QueryTruthTable(v -> v[1] & v[2], [
        (:isaromatic,), (:symbol, :N)
    ])
    and2 = QueryTruthTable(v -> v[1] & v[2] & v[3], [
        (:charge, 1), (:isaromatic,), (:symbol, :N)
    ])
    @test issubset(and1, symN)
    @test !issubset(symN, and1)
    @test issubset(and2, and1)
    @test !issubset(and1, and2)

    or1 = QueryTruthTable(v -> v[1] | v[2], [
        (:charge, 0), (:charge, 1)
    ])
    or2 = QueryTruthTable(v -> v[1] | v[2] | v[3], [
        (:charge, -1), (:charge, 0), (:charge, 1)
    ])
    @test issubset(ch1, or1)
    @test !issubset(or1, ch1)
    @test issubset(or1, or2)
    @test !issubset(or2, or1)

    nested1 = QueryTruthTable(
        v -> ~v[1] & (v[2] | v[3]) & (v[4] | v[5]),
        [
            (:isaromatic,), (:symbol, :O), (:symbol, :N),
            (:charge, 0), (:charge, 1)
        ]
    )
    nested2 = QueryTruthTable(
        v -> ~v[1] & (v[2] | v[3] | v[4]) & (v[5] | v[6] | v[7]),
        [
            (:isaromatic,), (:symbol, :O), (:symbol, :N), (:symbol, :S),
            (:charge, 0), (:charge, 1), (:smallest_ring, 6)
        ]
    )
    nested3 = QueryTruthTable(v -> ~v[1] & v[2] & (v[3] | v[4]), [
        (:isaromatic,), (:symbol, :O), (:charge, 0), (:charge, 1)
    ])
    @test issubset(nested1, nested2)
    @test !issubset(nested2, nested1)
    @test issubset(nested3, nested1)
    @test !issubset(nested1, nested3)

    # TODO: tautology
    # TODO: wildcard

    # TODO: disjoint varaibles (e.g. (A=2)|(A=3) is a subset of ~(A=1))
    """
    OorS = QueryTruthTable(v -> v[1] | v[2], [
        (:symbol, :O), (:symbol, :S)
    ])
    @test issubset(OorS, notN)
    @test !issubset(notN, OorS)
    """
    # TODO: recursive

end

"""
@testset "equivalence" begin

    rec1 = QueryFormula(:recursive, "[NH2]C")
    rec2 = QueryFormula(:recursive, "[NH2]CC")
    rec3 = QueryFormula(:recursive, "[NH2+0]C")
    rec4 = QueryFormula(:recursive, "[N,O;H2]C")
    and3 = QueryFormula(:and, Set([
        QueryFormula(:symbol, :N)
        QueryFormula(:or, Set([
            QueryFormula(:total_hydrogens, 2),
            QueryFormula(:total_hydrogens, 3)
        ]))
    ]))
    @test issubset(rec2, rec1)
    @test !issubset(rec1, rec2)
    @test issubset(rec3, rec1)
    @test !issubset(rec1, rec3)
    @test !issubset(rec4, rec1)
    @test issubset(rec1, rec4)
    @test issubset(rec1, and3)
    @test !issubset(rec4, and3)

    rec5 = QueryFormula(:recursive, "C-Cl")
    rec6 = QueryFormula(:or, Set([
        QueryFormula(:recursive, "C-Cl")
        QueryFormula(:recursive, "C-Br")
    ]))
    not4 = QueryFormula(:not, QueryFormula(:symbol, :O))
    @test issubset(rec5, rec6)
    @test !issubset(rec6, rec5)
    @test issubset(rec5, not4)
    @test !issubset(not4, rec5)


    fml4 = QueryFormula(:symbol, :S)
    fml5 = QueryFormula(:symbol, :N)
    and4 = QueryFormula(:and, Set([
        QueryFormula(:not, QueryFormula(:symbol, :C)),
        QueryFormula(:not, QueryFormula(:symbol, :N)),
        QueryFormula(:not, QueryFormula(:symbol, :O))
    ]))
    @test issubset(fml4, and4)
    @test !issubset(fml5, and4)

    or5 = QueryFormula(:or, Set([
        QueryFormula(:symbol, :N),
        QueryFormula(:symbol, :P)
    ]))
    and5 = QueryFormula(:and, Set([
        QueryFormula(:not, QueryFormula(:symbol, :C)),
        QueryFormula(:not, QueryFormula(:symbol, :S)),
        QueryFormula(:not, QueryFormula(:symbol, :O))
    ]))
    and6 = QueryFormula(:and, Set([
        QueryFormula(:symbol, :N),
        QueryFormula(:charge, 0),
        QueryFormula(:isaromatic, false)
    ]))
    and7 = QueryFormula(:and, Set([
        QueryFormula(:symbol, :C),
        QueryFormula(:charge, 0),
        QueryFormula(:isaromatic, false)
    ]))
    and8 = QueryFormula(:and, Set([
        QueryFormula(:or, Set([
            QueryFormula(:symbol, :C),
            QueryFormula(:symbol, :N)
        ])),
        QueryFormula(:isaromatic, false)
    ]))
    @test issubset(or5, and5)
    @test !issubset(and5, or5)
    @test issubset(and6, and5)
    @test !issubset(and7, and5)
    @test !issubset(and8, and5)

    and9 = QueryFormula(:and, Set([
        QueryFormula(:symbol, :F),
        QueryFormula(:isaromatic, false)
    ]))
    or6 = QueryFormula(:or, Set([
        QueryFormula(:symbol, :F),
        QueryFormula(:symbol, :Cl)
    ]))
    @test issubset(and9, or6)
    @test !issubset(or6, and9)
    and10 = QueryFormula(:and, Set([
        QueryFormula(:symbol, :C),
        QueryFormula(:isaromatic, false)
    ]))
    or7 = QueryFormula(:or, Set([
        QueryFormula(:and, Set([
            QueryFormula(:symbol, :C),
            QueryFormula(:isaromatic, false)
        ])),
        QueryFormula(:isaromatic, true)
    ]))
    @test issubset(and10, or7)
    @test !issubset(or7, and10)
end
"""

end # model.query