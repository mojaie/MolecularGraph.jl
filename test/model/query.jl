
@testset "model.query" begin

@testset "equivalence" begin
    @test QueryFormula(:v, :a) == QueryFormula(:v, :a)
    @test QueryFormula(:v, :b) != QueryFormula(:v, :c)
    @test QueryFormula(:x, 1) != QueryFormula(:v, 1)

    fml1 = QueryFormula(:x, "hoge")
    and1 = QueryFormula(:and, Set([
        QueryFormula(:x, "hoge"),
        QueryFormula(:y, "fuga"),
        QueryFormula(:z, "piyo")
    ]))
    and2 = QueryFormula(:and, Set([
        QueryFormula(:z, "piyo"),
        QueryFormula(:y, "fuga")
    ]))
    @test issubset(and1, fml1)
    @test !issubset(fml1, and1)
    @test and1 == and1
    @test and1 != and2
    @test issubset(and1, and2)  # seems not intuitive, but and1 result records should be less than and2
    @test !issubset(and2, and1)
    @test issubset(and1, QueryFormula(:x, "hoge"))
    @test !issubset(QueryFormula(:x, "hoge"), and1)

    fml2 = QueryFormula(:x, "fuga")
    or1 = QueryFormula(:or, Set([
        QueryFormula(:x, "hoge"),
        QueryFormula(:x, "fuga"),
        QueryFormula(:x, "piyo")
    ]))
    or2 = QueryFormula(:or, Set([
        QueryFormula(:x, "piyo"),
        QueryFormula(:x, "fuga")
    ]))
    @test issubset(fml2, or1)
    @test !issubset(or1, fml2)
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

    nota = QueryFormula(:not, QueryFormula(:v, :a))
    notb = QueryFormula(:not, QueryFormula(:v, :b))
    @test nota == nota
    @test nota != notb

    fml3 = QueryFormula(:x, 1)
    or3 = QueryFormula(:or, Set([
        QueryFormula(:x, 1),
        QueryFormula(:x, 2)
    ]))
    or4 = QueryFormula(:or, Set([
        QueryFormula(:v, 1),
        QueryFormula(:x, 1),
        QueryFormula(:x, 2)
    ]))
    not3 = QueryFormula(:not, QueryFormula(:x, 3))
    @test issubset(fml3, not3)
    @test !issubset(not3, fml3)
    @test or3 != not3
    @test issubset(or3, not3)
    @test !issubset(not3, or3)
    @test !issubset(or4, not3)

    rec1 = QueryFormula(:recursive, "[NH2]C")
    rec2 = QueryFormula(:recursive, "[NH2]CC")
    rec3 = QueryFormula(:recursive, "[NH2+0]C")
    rec4 = QueryFormula(:recursive, "[N,O;H2]C")
    and3 = QueryFormula(:and, Set([
        QueryFormula(:atomsymbol, :N)
        QueryFormula(:or, Set([
            QueryFormula(:hydrogenconnected, 2),
            QueryFormula(:hydrogenconnected, 3)
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
    not4 = QueryFormula(:not, QueryFormula(:atomsymbol, :O))
    @test issubset(rec5, rec6)
    @test !issubset(rec6, rec5)
    @test issubset(rec5, not4)
    @test !issubset(not4, rec5)


    fml4 = QueryFormula(:atomsymbol, :S)
    fml5 = QueryFormula(:atomsymbol, :N)
    and4 = QueryFormula(:and, Set([
        QueryFormula(:not, QueryFormula(:atomsymbol, :C)),
        QueryFormula(:not, QueryFormula(:atomsymbol, :N)),
        QueryFormula(:not, QueryFormula(:atomsymbol, :O))
    ]))
    @test issubset(fml4, and4)
    @test !issubset(fml5, and4)

    or5 = QueryFormula(:or, Set([
        QueryFormula(:atomsymbol, :N),
        QueryFormula(:atomsymbol, :P)
    ]))
    and5 = QueryFormula(:and, Set([
        QueryFormula(:not, QueryFormula(:atomsymbol, :C)),
        QueryFormula(:not, QueryFormula(:atomsymbol, :S)),
        QueryFormula(:not, QueryFormula(:atomsymbol, :O))
    ]))
    and6 = QueryFormula(:and, Set([
        QueryFormula(:atomsymbol, :N),
        QueryFormula(:charge, 0),
        QueryFormula(:isaromatic, false)
    ]))
    and7 = QueryFormula(:and, Set([
        QueryFormula(:atomsymbol, :C),
        QueryFormula(:charge, 0),
        QueryFormula(:isaromatic, false)
    ]))
    and8 = QueryFormula(:and, Set([
        QueryFormula(:or, Set([
            QueryFormula(:atomsymbol, :C),
            QueryFormula(:atomsymbol, :N)
        ])),
        QueryFormula(:isaromatic, false)
    ]))
    @test issubset(or5, and5)
    @test !issubset(and5, or5)
    @test issubset(and6, and5)
    @test !issubset(and7, and5)
    @test !issubset(and8, and5)

    and9 = QueryFormula(:and, Set([
        QueryFormula(:atomsymbol, :F),
        QueryFormula(:isaromatic, false)
    ]))
    or6 = QueryFormula(:or, Set([
        QueryFormula(:atomsymbol, :F),
        QueryFormula(:atomsymbol, :Cl)
    ]))
    @test issubset(and9, or6)
    @test !issubset(or6, and9)
    and10 = QueryFormula(:and, Set([
        QueryFormula(:atomsymbol, :C),
        QueryFormula(:isaromatic, false)
    ]))
    or7 = QueryFormula(:or, Set([
        QueryFormula(:and, Set([
            QueryFormula(:atomsymbol, :C),
            QueryFormula(:isaromatic, false)
        ])),
        QueryFormula(:isaromatic, true)
    ]))
    @test issubset(and10, or7)
    @test !issubset(or7, and10)
end

@testset "tidyformula" begin
    and1 = QueryFormula(:and, Set([
        QueryFormula(:x, :a),
        QueryFormula(:and, Set([
            QueryFormula(:y, :b),
            QueryFormula(:z, :c)
        ]))
    ]))
    @test tidyformula(and1) == QueryFormula(:and, Set([
        QueryFormula(:x, :a),
        QueryFormula(:y, :b),
        QueryFormula(:z, :c)
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
            QueryFormula(:x, :a),
            QueryFormula(:y, :b)
        ])),
        QueryFormula(:and, Set([
            QueryFormula(:x, :c),
            QueryFormula(:y, :d)
        ])),
        QueryFormula(:or, Set([
            QueryFormula(:x, :e),
            QueryFormula(:y, :f)
        ]))
    ]))
    @test tidyformula(or2) == QueryFormula(:or, Set([
        QueryFormula(:and, Set([
            QueryFormula(:x, :a),
            QueryFormula(:y, :b)
        ])),
        QueryFormula(:and, Set([
            QueryFormula(:x, :c),
            QueryFormula(:y, :d)
        ])),
        QueryFormula(:x, :e),
        QueryFormula(:y, :f)
    ]))

    not1 = QueryFormula(:not, QueryFormula(:v, :a))
    @test tidyformula(not1) == not1

    and3 = QueryFormula(:and, Set([
        QueryFormula(:and, Set([
            QueryFormula(:and, Set([
                QueryFormula(:a, 1),
                QueryFormula(:and, Set([
                    QueryFormula(:b, 1),
                    QueryFormula(:c, 1)
                ]))
            ])),
            QueryFormula(:d, 1)
        ])),
        QueryFormula(:e, 1)
    ]))
    @test tidyformula(and3) == QueryFormula(:and, Set([
        QueryFormula(:a, 1),
        QueryFormula(:b, 1),
        QueryFormula(:c, 1),
        QueryFormula(:d, 1),
        QueryFormula(:e, 1)
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

    dist3 = QueryFormula(:and, Set([
        QueryFormula(:hydrogenconnected, 1),
        QueryFormula(:or, Set([
            QueryFormula(:and, Set([
                QueryFormula(:atomsymbol, :O),
                QueryFormula(:isaromatic, false)
            ])),
            QueryFormula(:and, Set([
                QueryFormula(:atomsymbol, :N),
                QueryFormula(:isaromatic, false)
            ]))
        ]))
    ]))
    @test tidyformula(dist3) == QueryFormula(:and, Set([
        QueryFormula(:hydrogenconnected, 1),
        QueryFormula(:or, Set([
            QueryFormula(:atomsymbol, :O),
            QueryFormula(:atomsymbol, :N)
        ])),
        QueryFormula(:isaromatic, false)
    ]))

    abs3 = QueryFormula(:and, Set([
        QueryFormula(:any, true),
        QueryFormula(:a, 1)
    ]))
    @test tidyformula(abs3) == QueryFormula(:a, 1)

    abs4 = QueryFormula(:or, Set([
        QueryFormula(:any, true),
        QueryFormula(:a, 1)
    ]))
    @test tidyformula(abs4) == QueryFormula(:any, true)

    abs5 = QueryFormula(:and, Set([
        QueryFormula(:any, false),
        QueryFormula(:a, 1)
    ]))
    @test tidyformula(abs5) == QueryFormula(:any, false)

    abs6 = QueryFormula(:or, Set([
        QueryFormula(:any, false),
        QueryFormula(:a, 1)
    ]))
    @test tidyformula(abs6) == QueryFormula(:a, 1)

    abs7 = QueryFormula(:and, Set([
        QueryFormula(:or, Set([
            QueryFormula(:or, Set([
                QueryFormula(:any, true),
                QueryFormula(:c, 3)
            ])),
            QueryFormula(:b, 2)
        ])),
        QueryFormula(:a, 1)
    ]))
    @test tidyformula(abs7) == QueryFormula(:a, 1)

    disjoint1 = QueryFormula(:and, Set([
        QueryFormula(:v, :a)
        QueryFormula(:v, :c)
        QueryFormula(:v, :e)
    ]))
    @test tidyformula(disjoint1) == QueryFormula(:any, false)

    disjoint2 = QueryFormula(:and, Set([
        QueryFormula(:not, QueryFormula(:v, :a))
        QueryFormula(:not, QueryFormula(:v, :c))
        QueryFormula(:not, QueryFormula(:v, :e))
    ]))
    @test tidyformula(disjoint2) == disjoint2

    disjoint3 = QueryFormula(:and, Set([
        QueryFormula(:or, Set([
            QueryFormula(:x, 1),
            QueryFormula(:x, 2)
        ])),
        QueryFormula(:or, Set([
            QueryFormula(:x, 2),
            QueryFormula(:x, 3)
        ])),
        QueryFormula(:or, Set([
            QueryFormula(:x, 3),
            QueryFormula(:x, 1)
        ])),
    ]))
    @test tidyformula(disjoint3) == QueryFormula(:any, false)

    disjoint4 = QueryFormula(:and, Set([
        QueryFormula(:or, Set([
            QueryFormula(:x, 1),
            QueryFormula(:y, :a)
        ])),
        QueryFormula(:or, Set([
            QueryFormula(:x, 2),
            QueryFormula(:y, :b)
        ])),
        QueryFormula(:or, Set([
            QueryFormula(:x, 3),
            QueryFormula(:y, :c)
        ])),
    ]))
    @test tidyformula(disjoint4) == disjoint4
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