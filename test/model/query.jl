
using MolecularGraph:
    querypropmap, generate_queryfunc, querymatch, optimize_query

@testset "model.query" begin

@testset "truthtable" begin
    # tautology
    """
    q \\ r  T other F
    T      T   F   F
    other  T   ?   F
    F      T   T   T
    """
    props = [(:symbol, :C)]
    t = QueryTruthTable(v -> true, props)
    f = QueryTruthTable(v -> false, props)
    o = QueryTruthTable(v -> v[1], props)
    @test t == t
    @test f == f
    @test !issubset(t, o)
    @test !issubset(t, f)
    @test issubset(o, t)
    @test !issubset(o, f)
    @test issubset(f, t)
    @test issubset(f, o)

    props = [(:symbol, :C), (:symbol, :N)]
    symC = QueryTruthTable(v -> v[1], props)
    symN = QueryTruthTable(v -> v[2], props)
    notN = QueryTruthTable(v -> ~v[2], props)
    props = [(:total_hydrogens, 1), (:charge, 1)]
    hc1 = QueryTruthTable(v -> v[1], props)
    ch1 = QueryTruthTable(v -> v[2], props)
    @test symC == symC
    @test symC != symN
    @test symN != symC
    @test symN != notN
    @test notN != symN
    @test hc1 == hc1
    @test hc1 != ch1
    @test ch1 != hc1

    props = [(:charge, 1), (:isaromatic,), (:symbol, :N)]
    and1 = QueryTruthTable(v -> v[2] & v[3], props)
    and2 = QueryTruthTable(v -> v[1] & v[2] & v[3], props)
    @test issubset(and2, and1)
    @test !issubset(and1, and2)

    props = [(:charge, -1), (:charge, 0), (:charge, 1)]
    or1 = QueryTruthTable(v -> v[2] | v[3], props)
    or2 = QueryTruthTable(v -> v[1] | v[2] | v[3], props)
    @test issubset(or1, or2)
    @test !issubset(or2, or1)

    props = [
        (:charge, 0), (:charge, 1), (:isaromatic,),
        (:smallest_ring, 6), (:symbol, :N), (:symbol, :O), (:symbol, :S)
    ]
    nested1 = QueryTruthTable(
        v -> ~v[3] & (v[5] | v[6]) & (v[1] | v[2]), props)
    nested2 = QueryTruthTable(
        v -> ~v[3] & (v[5] | v[6] | v[7]) & (v[1] | v[2] | v[4]), props)
    nested3 = QueryTruthTable(
        v -> ~v[3] & v[6] & (v[1] | v[2]), props)
    @test issubset(nested1, nested2)
    @test !issubset(nested2, nested1)
    @test issubset(nested3, nested1)
    @test !issubset(nested1, nested3)
end

@testset "querynode" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    @test QueryLiteral(:a, 2) > QueryLiteral(:a, 1)
    @test QueryLiteral(:a, 2) < QueryLiteral(:b, 2)
    @test ([QueryLiteral(:a, 1), QueryLiteral(:b, 2), QueryLiteral(:c, 3)]
            == [QueryLiteral(:a, 1), QueryLiteral(:b, 2), QueryLiteral(:c, 3)])

    isc = generate_queryfunc(
        QueryLiteral(:symbol, :C), [QueryLiteral(:symbol, :C)])
    @test isc([true])
    notc = generate_queryfunc(
        QueryOperator(:not, [QueryLiteral(:symbol, :C)]), [QueryLiteral(:symbol, :C)])
    @test notc([false])
    af = generate_queryfunc(QueryAny(false), [QueryLiteral(:symbol, :C)])
    @test !af([true])
    @test !af([false])

    op1 = QueryOperator(:and, [
        QueryLiteral(:symbol, :C),
        QueryOperator(:not, [QueryLiteral(:isaromatic)]),
        QueryOperator(:or, [
            QueryLiteral(:connectivity, 3),
            QueryLiteral(:connectivity, 4)
        ])
    ])
    @test length(querypropmap(op1)) == 3
    tbl1 = QueryTruthTable(op1)
    func1 = QueryTruthTable(
        v -> v[4] & ~v[3] & (v[1] | v[2]),
        [(:connectivity, 3), (:connectivity, 4), (:isaromatic,), (:symbol, :C)]
    )
    @test querymatch(tbl1, func1, true)
    # global_logger(default_logger)
end

@testset "optimize" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    c1 = QueryOperator(:and, [QueryLiteral(:symbol, :C), QueryAny(true)])
    @test optimize_query(c1) == QueryLiteral(:symbol, :C)
    c2 = QueryOperator(:or, [QueryLiteral(:symbol, :C), QueryAny(true)])
    @test optimize_query(c2) == QueryAny(true)
    c3 = QueryOperator(:and, [QueryAny(false), QueryLiteral(:symbol, :C)])
    @test optimize_query(c3) == QueryAny(false)
    c4 = QueryOperator(:or, [QueryAny(false), QueryLiteral(:symbol, :C)])
    @test optimize_query(c4) == QueryLiteral(:symbol, :C)

    notn = QueryOperator(:not, [
        QueryOperator(:and, [
            QueryLiteral(:symbol, :N),
            QueryOperator(:not, [QueryLiteral(:isaromatic)])
        ])
    ])
    notn_ = optimize_query(notn)
    @test notn_.key === :or
    @test notn_.value[1].key === :not
    @test notn_.value[2] == QueryLiteral(:isaromatic)
    # global_logger(default_logger)
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