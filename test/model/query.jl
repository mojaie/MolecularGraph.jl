
@testset "model.query" begin

@testset "equivalence" begin
    @test QueryFormula(:v, :a) == QueryFormula(:v, :a)
    @test QueryFormula(:v, :b) != QueryFormula(:v, :c)
    @test QueryFormula(:x, 1) != QueryFormula(:v, 1)

    nota = QueryFormula(:not, QueryFormula(:v, :a))
    notb = QueryFormula(:not, QueryFormula(:v, :b))
    @test nota == nota
    @test nota != notb

    and1 = QueryFormula(:and, Set([
        QueryFormula(:x, "hoge"),
        QueryFormula(:y, "fuga"),
        QueryFormula(:z, "piyo")
    ]))
    and2 = QueryFormula(:and, Set([
        QueryFormula(:z, "piyo"),
        QueryFormula(:y, "fuga")
    ]))
    @test and1 == and1
    @test and1 != and2
    @test issubset(and1, and2)  # seems not intuitive, but and1 result records should be less than and2
    @test !issubset(and2, and1)
    @test issubset(and1, QueryFormula(:x, "hoge"))
    @test !issubset(QueryFormula(:x, "hoge"), and1)

    or1 = QueryFormula(:or, Set([
        QueryFormula(:x, "hoge"),
        QueryFormula(:x, "fuga"),
        QueryFormula(:x, "piyo")
    ]))
    or2 = QueryFormula(:or, Set([
        QueryFormula(:x, "piyo"),
        QueryFormula(:x, "fuga")
    ]))
    @test or1 == or1
    @test or1 != or2
    @test issubset(or2, or1)  # opposit to the `and` case
    @test !issubset(or1, or2)

    nested1 = QueryFormula(:and, Set([
        QueryFormula(:x, 1),
        QueryFormula(:or, Set([
            QueryFormula(:v, :a),
            QueryFormula(:v, :b)
        ])),
        QueryFormula(:or, Set([
            QueryFormula(:n, true),
            QueryFormula(:m, true)
        ]))
    ]))
    nested2 = QueryFormula(:and, Set([
        QueryFormula(:x, 1),
        QueryFormula(:or, Set([
            QueryFormula(:v, :a),
            QueryFormula(:v, :b),
            QueryFormula(:v, :c)
        ])),
        QueryFormula(:or, Set([
            QueryFormula(:n, true),
            QueryFormula(:m, true),
            QueryFormula(:o, true)
        ]))
    ]))
    nested3 = QueryFormula(:and, Set([
        QueryFormula(:x, 1),
        QueryFormula(:v, :a),
        QueryFormula(:or, Set([
            QueryFormula(:n, true),
            QueryFormula(:m, true)
        ]))
    ]))
    @test nested1 == nested1
    @test nested1 != nested2
    @test issubset(nested1, nested2)
    @test !issubset(nested2, nested1)
    @test issubset(nested3, nested1)
    @test !issubset(nested1, nested3)
end

@testset "tidyformula" begin
    and1 = QueryFormula(:and, Set([
        QueryFormula(:v, :a),
        QueryFormula(:and, Set([
            QueryFormula(:v, :b),
            QueryFormula(:v, :c)
        ]))
    ]))
    @test tidyformula(and1) == QueryFormula(:and, Set([
        QueryFormula(:v, :a),
        QueryFormula(:v, :b),
        QueryFormula(:v, :c)
    ]))

    or1 = QueryFormula(:or, Set([
        QueryFormula(:or, Set([
            QueryFormula(:x, :a),
            QueryFormula(:y, :a)
        ])),
        QueryFormula(:z, :a)
    ]))
    @test tidyformula(or1) == QueryFormula(:or, Set([
        QueryFormula(:x, :a),
        QueryFormula(:y, :a),
        QueryFormula(:z, :a)
    ]))

    and2 = QueryFormula(:and, Set([
        QueryFormula(:and, Set([
            QueryFormula(:x, 1),
            QueryFormula(:y, 1)
        ])),
        QueryFormula(:v, :a),
        QueryFormula(:or, Set([
            QueryFormula(:z, 0),
            QueryFormula(:z, nothing)
        ])),
        QueryFormula(:not, QueryFormula(:a, "hoge"))
    ]))
    @test tidyformula(and2) == QueryFormula(:and, Set([
        QueryFormula(:v, :a),
        QueryFormula(:x, 1),
        QueryFormula(:y, 1),
        QueryFormula(:or, Set([
            QueryFormula(:z, 0),
            QueryFormula(:z, nothing)
        ])),
        QueryFormula(:not, QueryFormula(:a, "hoge"))
    ]))

    or2 = QueryFormula(:or, Set([
        QueryFormula(:and, Set([
            QueryFormula(:v, :a),
            QueryFormula(:v, :b)
        ])),
        QueryFormula(:and, Set([
            QueryFormula(:v, :c),
            QueryFormula(:v, :d)
        ])),
        QueryFormula(:or, Set([
            QueryFormula(:v, :e),
            QueryFormula(:v, :f)
        ]))
    ]))
    @test tidyformula(or2) == QueryFormula(:or, Set([
        QueryFormula(:and, Set([
            QueryFormula(:v, :a),
            QueryFormula(:v, :b)
        ])),
        QueryFormula(:and, Set([
            QueryFormula(:v, :c),
            QueryFormula(:v, :d)
        ])),
        QueryFormula(:v, :e),
        QueryFormula(:v, :f)
    ]))

    not1 = QueryFormula(:not, QueryFormula(:v, :a))
    @test tidyformula(not1) == not1

    and3 = QueryFormula(:and, Set([
        QueryFormula(:and, Set([
            QueryFormula(:and, Set([
                QueryFormula(:v, :a),
                QueryFormula(:and, Set([
                    QueryFormula(:v, :b),
                    QueryFormula(:v, :c)
                ]))
            ])),
            QueryFormula(:v, :d)
        ])),
        QueryFormula(:v, :e)
    ]))
    @test tidyformula(and3) == QueryFormula(:and, Set([
        QueryFormula(:v, :a),
        QueryFormula(:v, :b),
        QueryFormula(:v, :c),
        QueryFormula(:v, :d),
        QueryFormula(:v, :e)
    ]))

    abs1 = QueryFormula(:and, Set([
        QueryFormula(:x, "important"),
        QueryFormula(:or, Set([
            QueryFormula(:x, "important"),
            QueryFormula(:z, "optional")
        ]))
    ]))
    @test tidyformula(abs1) == QueryFormula(:x, "important")

    abs2 = QueryFormula(:and, Set([
        QueryFormula(:or, Set([
            QueryFormula(:a, 1),
            QueryFormula(:a, 2),
            QueryFormula(:a, 3)
        ]))
        QueryFormula(:or, Set([
            QueryFormula(:a, 1),
            QueryFormula(:a, 2)
        ]))
    ]))
    @test tidyformula(abs2) == QueryFormula(:or, Set([
        QueryFormula(:a, 1),
        QueryFormula(:a, 2)
    ]))

    dist1 = QueryFormula(:or, Set([
        QueryFormula(:and, Set([
            QueryFormula(:a, 1),
            QueryFormula(:b, 2),
            QueryFormula(:c, 3),
        ])),
        QueryFormula(:and, Set([
            QueryFormula(:a, 1),
            QueryFormula(:b, 2),
            QueryFormula(:c, 4)
        ]))
    ]))
    @test tidyformula(dist1) == QueryFormula(:and, Set([
        QueryFormula(:a, 1),
        QueryFormula(:b, 2),
        QueryFormula(:or, Set([
            QueryFormula(:c, 3)
            QueryFormula(:c, 4)
        ]))
    ]))

    dist2 = QueryFormula(:and, Set([
        QueryFormula(:or, Set([
            QueryFormula(:x, false),
            QueryFormula(:y, nothing)
        ])),
        QueryFormula(:or, Set([
            QueryFormula(:x, false),
            QueryFormula(:z, 3),
            QueryFormula(:z, 4)
        ]))
    ]))
    @test tidyformula(dist2) == QueryFormula(:or, Set([
        QueryFormula(:x, false),
        QueryFormula(:and, Set([
            QueryFormula(:y, nothing)
            QueryFormula(:or, Set([
                QueryFormula(:z, 3),
                QueryFormula(:z, 4)
            ]))
        ]))
    ]))
end

@testset "querymol" begin
    q = querymol(SmartsAtom, SmartsBond)
    addatom!(q, SmartsAtom(QueryFormula(:any, true)))
    addatom!(q, SmartsAtom(QueryFormula(:any, true)))
    addbond!(q, 1, 2, SmartsBond(QueryFormula(:any, true)))
    @test atomcount(q) == 2
    @test bondcount(q) == 1
    @test hasbond(q, 1, 2)
end

end # model.query