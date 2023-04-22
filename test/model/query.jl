
using MolecularGraph:
    querypropmap, generate_queryfunc, querymatch, optimize_query,
    fragment!,specialize_nonaromatic!, remove_hydrogens!

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
    @test QueryAny(true) == QueryAny(true)
    @test QueryAny(true) != QueryAny(false)
    @test QueryLiteral(:a, 1) == QueryLiteral(:a, 1)
    @test QueryLiteral(:a, "hoge") != QueryLiteral(:a, "fuga")
    @test QueryLiteral(:charge, 2) != QueryLiteral(:total_hydrogens, 2)
    @test QueryOperator(:and, [
        QueryLiteral(:a, 1),
        QueryOperator(:not, [QueryLiteral(:a, 1)]),
        QueryOperator(:or, [QueryLiteral(:b, :hoge), QueryLiteral(:c, :fuga)]),

    ]) == QueryOperator(:and, [
        QueryOperator(:not, [QueryLiteral(:a, 1)]),
        QueryOperator(:or, [QueryLiteral(:c, :fuga), QueryLiteral(:b, :hoge)]),
        QueryLiteral(:a, 1)
    ])
    @test QueryOperator(:or, [
        QueryOperator(:or, [
            QueryOperator(:or, [
                QueryOperator(:or, [QueryLiteral(:a, true), QueryLiteral(:b, false)]),
                QueryLiteral(:c, false)
            ]),
            QueryLiteral(:d, false)
        ])
    ]) != QueryOperator(:or, [
        QueryOperator(:or, [
            QueryOperator(:or, [
                QueryOperator(:or, [QueryLiteral(:a, true), QueryLiteral(:b, true)]),
                QueryLiteral(:c, false)
            ]),
            QueryLiteral(:d, false)
        ])
    ])

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
    @test optimize_query(notn) == QueryOperator(:or, [
        QueryOperator(:not, [QueryLiteral(:symbol, :N)]),
        QueryLiteral(:isaromatic)
    ])
    # global_logger(default_logger)
end

@testset "pains" begin
    narom = QueryOperator(:not, [QueryLiteral(:isaromatic)])
    state = SMARTSParser{MolGraph{Int,QueryTree,QueryTree}}(
        "n1(-[#6])c(c(-[#1])c(c1-[#6]=[#7]-[#7])-[#1])-[#1]")  # hzone_pyrrol(64)
    fragment!(state)
    pains1 = MolGraph{Int,QueryTree,QueryTree}(
        state.edges, state.vprops, state.eprops, gprop_map=Dict(:connectivity => state.connectivity))
    specialize_nonaromatic!(pains1)
    @test get_prop(pains1, 2, :tree) == QueryLiteral(:symbol, :C)  # -[#6] still can be aromatic
    @test get_prop(pains1, 5, :tree) == QueryLiteral(:symbol, :H)  # non-aromatic symbols are not affected
    @test get_prop(pains1, 8, :tree) == QueryOperator(:and, [QueryLiteral(:symbol, :C), narom])
    @test get_prop(pains1, 9, :tree) == QueryOperator(:and, [QueryLiteral(:symbol, :N), narom])
    @test get_prop(pains1, 10, :tree) == QueryOperator(:and, [QueryLiteral(:symbol, :N), narom])
    remove_hydrogens!(pains1)
    @test nv(pains1) == 9
    @test degree(pains1, 3) == 2
    @test degree(pains1, 4) == 2
    @test degree(pains1, 6) == 2

    pains2 = smartstomol("[!#6&!#1]=[#6]1[#6]=,:[#6][#6](=[!#6&!#1])[#6]=,:[#6]1")  # quinone_A(370)
    @test get_prop(pains2, 1, :tree) == QueryOperator(:not, [QueryLiteral(:symbol, :C)])
    @test get_prop(pains2, 2, :tree) == QueryOperator(:and, [QueryLiteral(:symbol, :C), narom])
    @test get_prop(pains2, 4, :tree) == QueryLiteral(:symbol, :C)
    @test get_prop(pains2, 6, :tree) == QueryOperator(:not, [QueryLiteral(:symbol, :C)])
end

end # model.query