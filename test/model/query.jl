
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

end # model.query