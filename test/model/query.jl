
@testset "model.query" begin

@testset "querynode" begin
    @test qand() == qand()
    @test qand() != qor()
    @test qnot() == qnot()
    @test qanytrue() != qnot()
    @test qtrue(:x) == qtrue(:x)
    @test qtrue(:x) != qtrue(:y)
    @test qeq(:x, "1") == qeq(:x, "1")
    @test qeq(:x, "1") != qeq(:x, "2")
    @test qeq(:x, "1") != qeq(:y, "1")
    @test hash(qor(), hash(qanytrue())) == hash(qor(), hash(qanytrue()))
    @test hash(qand(), hash(qnot())) != hash(qand(), hash(qand()))
    @test hash(qtrue(:x), hash(qeq(:y, "true"))) == hash(qtrue(:x), hash(qeq(:y, "true")))
    @test hash(qtrue(:x), hash(qeq(:y, "true"))) != hash(qtrue(:x), hash(qeq(:y, "false")))
end

@testset "querytree" begin

    a = QueryAtom(
        [(4, 2), (4, 3), (4, 1), (3, 5), (1, 6), (1, 7)],
        [qand(), qeq(:a, "1"), qnot(), qor(), qeq(:a, "1"), qeq(:b, "hoge"), qeq(:c, "2")]
    )
    @test canonical(a)[2] == [4, 1, 6, 7, 2, 3, 5]

    @test QueryAtom() == QueryAtom()
    @test QueryAtom() != QueryBond()
    @test hash(QueryAtom()) == hash(QueryAtom())

    @test QueryAtom(
        [(1, 2), (1, 3), (1, 4)],
        [qand(), qeq(:b, "1"), qeq(:a, "3"), qeq(:b, "2")]
    ) == QueryAtom(
        [(4, 1), (4, 2), (4, 3)],
        [qeq(:b, "1"), qeq(:a, "3"), qeq(:b, "2"), qand()]
    )

    @test QueryBond(
        [(1, 2), (1, 3), (1, 4), (3, 5), (4, 6), (4, 7)],
        [qor(), qeq(:a, "1"), qnot(), qand(), qeq(:a, "1"), qeq(:b, "hoge"), qeq(:c, "fuga")]
    ) == QueryBond(
        [(1, 2), (1, 3), (1, 4), (3, 5), (4, 6), (4, 7)],
        [qor(), qeq(:a, "1"), qnot(), qand(), qeq(:a, "1"), qeq(:c, "fuga"), qeq(:b, "hoge")]
    )

    @test QueryAtom(
        [(1, 2), (1, 3), (1, 4), (3, 5), (4, 6), (4, 7)],
        [qor(), qeq(:a, "1"), qnot(), qand(), qeq(:a, "1"), qeq(:b, "hoge"), qeq(:c, "1")]
    ) != QueryAtom(
        [(1, 2), (1, 3), (1, 4), (3, 5), (4, 6), (4, 7)],
        [qor(), qeq(:a, "1"), qnot(), qand(), qeq(:a, "1"), qeq(:b, "hoge"), qeq(:c, "2")]
    )

    x = QueryAtom(
        [(1, 2), (1, 3), (1, 4), (2, 5), (2, 6), (3, 7), (3, 8), (4, 9), (4, 10)],
        [qand(), qor(), qor(), qor(), qeq(:c, "1"), qeq(:c, "2"),
        qeq(:b, "1"), qeq(:b, "2"), qeq(:a, "1"), qeq(:a, "2")]
    )
    y = QueryAtom(
        [(1, 2), (1, 3), (1, 4), (2, 5), (2, 6), (3, 8), (3, 7), (4, 9), (4, 10)],
        [qand(), qor(), qor(), qor(), qeq(:a, "1"), qeq(:a, "2"),
        qeq(:b, "1"), qeq(:b, "2"), qeq(:c, "1"), qeq(:c, "2")]
    )
    @test x == y
    @test hash(x) == hash(y)

    # edit
    qa = QueryAtom(
        [(1, 2), (1, 3), (1, 4)],
        [qand(), qeq(:b, "1"), qeq(:a, "3"), qeq(:b, "2")]
    )
    @test root(qa) == 1
    add_qnode!(qa, 1, qeq(:c, "3"))
    newn = add_qnode!(qa, qtrue(:x))
    newa = add_qnode!(qa, qand())
    add_qedge!(qa, newa, newn)
    add_qedge!(qa, newa, 1)
    @test qa == QueryAtom(
        [(1, 2), (1, 3), (1, 4), (1, 5), (7, 6), (7, 1)],
        [qand(), qeq(:b, "1"), qeq(:a, "3"), qeq(:b, "2"),
        qeq(:c, "3"), qtrue(:x), qand()]
    )
    @test root(qa) == 7
    rem_qnode!(qa, 3)
    rem_qnodes!(qa, [3, 6])
    set_qnode!(qa, 3, qtrue(:x))
    @test qa == QueryAtom(
        [(1, 2), (1, 3), (1, 4)],
        [qand(), qeq(:b, "1"), qtrue(:x), qeq(:b, "2")]
    )
end

@testset "evaluate" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    null = generate_queryfunc(QueryAtom(), QueryNode[])
    @test !null([])
    isc = generate_queryfunc(
        QueryAtom(Tuple{Int,Int}[], [qeq(:symbol, "C")]),
        [qeq(:symbol, "C")]
    )
    @test isc([true])
    notc = generate_queryfunc(
        QueryAtom([(1, 2)], [qnot(), qeq(:symbol, "C")]),
        [qeq(:symbol, "C")]
    )
    @test notc([false])
    af = generate_queryfunc(
        QueryAtom([(1, 2)], [qnot(), qanytrue()]),
        QueryNode[]
    )
    @test !af([true])
    @test !af([false])

    qa = generate_queryfunc(
        QueryAtom(
            [(1, 2), (1, 3), (3, 4)],
            [qand(), qeq(:symbol, "N"), qnot(), qtrue(:isaromatic)]),
        [qeq(:symbol, "N"), qeq(:x, "1"), qtrue(:isaromatic)]
    )
    @test !qa([true, true, true])
    @test !qa([true, false, true])
    @test !qa([false, false, true])
    @test qa([true, false, false])
    # global_logger(default_logger)
end

@testset "preprocess" begin
    # [#6]-[#6](-[#7])=[#6] -> [#6]-C(-N)=C
    atoms = [
        QueryAtom(Tuple{Int,Int}[], [qeq(:symbol, "C")]),
        QueryAtom(Tuple{Int,Int}[], [qeq(:symbol, "C")]),
        QueryAtom(Tuple{Int,Int}[], [qeq(:symbol, "N")]),
        QueryAtom(Tuple{Int,Int}[], [qeq(:symbol, "C")])
    ]
    bonds = [
        QueryBond([(1, 2), (1, 3), (3, 4)], [qand(),qeq(:order, "1"), qnot(), qtrue(:isaromatic)]),
        QueryBond([(1, 2), (1, 3), (3, 4)], [qand(),qeq(:order, "1"), qnot(), qtrue(:isaromatic)]),
        QueryBond([(1, 2), (1, 3), (3, 4)], [qand(),qeq(:order, "2"), qnot(), qtrue(:isaromatic)])
    ]
    qmol = MolGraph(Edge.([(1, 2), (2, 3), (2, 4)]), atoms, bonds)
    specialize_nonaromatic!(qmol)
    @test qmol.vprops[1] == QueryAtom(Tuple{Int,Int}[], [qeq(:symbol, "C")])
    @test qmol.vprops[3] == QueryAtom(
        [(1, 2), (1, 3), (3, 4)], [qand(), qeq(:symbol, "N"), qnot(), qtrue(:isaromatic)])
    @test qmol.vprops[4] == QueryAtom(
        [(1, 2), (1, 3), (3, 4)], [qand(), qeq(:symbol, "C"), qnot(), qtrue(:isaromatic)])

    noth = QueryAtom([(2, 1)], [qeq(:symbol, "H"), qnot()])
    resolve_not_hydrogen!(noth)
    @test noth == QueryAtom(Tuple{Int,Int}[], [qanytrue()])

    # [#1][!#1]([#1])[#1] -> [*]
    atoms = [
        QueryAtom(Tuple{Int,Int}[], [qeq(:symbol, "H")]),
        QueryAtom([(1, 2)], [qnot(), qeq(:symbol, "H")]),
        QueryAtom(Tuple{Int,Int}[], [qeq(:symbol, "H")]),
        QueryAtom(Tuple{Int,Int}[], [qeq(:symbol, "H")])
    ]
    bonds = [
        QueryBond(
            [(1, 2), (2, 3), (2, 4), (4, 5), (1, 6)],
            [qor(), qand(), qeq(:order, "1"), qnot(), qtrue(:isaromatic), qtrue(:isaromatic)]
        ) for _ in 1:3
    ]
    qmol = MolGraph(Edge.([(2, 1), (2, 3), (2, 4)]), atoms, bonds)
    remove_hydrogens!(qmol)
    @test nv(qmol) == 1
    @test qmol.vprops[1] == QueryAtom(
        [(1, 2), (1, 3), (1, 4), (1, 5), (3, 6), (4, 7), (5, 8)],
        [qand(), qanytrue(), qnot(), qnot(), qnot(), qeq(:total_hydrogens, "0"),
        qeq(:total_hydrogens, "1"), qeq(:total_hydrogens, "2")])
end


"""
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
"""

"""
@testset "pains" begin
    narom = QueryOperator(:not, [QueryLiteral(:isaromatic)])
    state = SMARTSParser{SMARTSMolGraph}(
        "n1(-[#6])c(c(-[#1])c(c1-[#6]=[#7]-[#7])-[#1])-[#1]")  # hzone_pyrrol(64)
    fragment!(state)
    pains1 = MolGraph{SMARTSMolGraph}(
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
"""

end # model.query