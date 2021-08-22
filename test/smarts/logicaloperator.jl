
using MolecularGraph:
    forward!, lgnot!, lghighand!, lgor!, lglowand!


@testset "smarts.logicaloperator" begin

@testset "logicaloperator" begin
    # Test function
    function parseabc(state)
        if read(state) == 'a'
            forward!(state)
            return QueryFormula(:v, :a)
        elseif read(state) == 'b'
            forward!(state)
            return QueryFormula(:v, :b)
        elseif read(state) == 'c'
            forward!(state)
            return QueryFormula(:v, :c)
        end
    end

    state = SmartsParser("!a", false)
    not1 = lgnot!(state, parseabc)
    @test not1 == QueryFormula(:not, QueryFormula(:v, :a))

    state = SmartsParser("a&b", false)
    and1 = lghighand!(state, parseabc)
    @test and1 == QueryFormula(:and, Set([
        QueryFormula(:v, :a),
        QueryFormula(:v, :b)
    ]))

    state = SmartsParser("abc", false)
    and2 = lghighand!(state, parseabc)
    @test and2 == QueryFormula(:and, Set([
        QueryFormula(:v, :a),
        QueryFormula(:v, :b),
        QueryFormula(:v, :c)
    ]))

    state = SmartsParser("!ab", false)
    not2 = lghighand!(state, parseabc)
    @test not2 == QueryFormula(:and, Set([
        QueryFormula(:not, QueryFormula(:v, :a)),
        QueryFormula(:v, :b)
    ]))

    state = SmartsParser("a,b", false)
    or1 = lgor!(state, parseabc)
    @test or1 == QueryFormula(:or, Set([
        QueryFormula(:v, :a),
        QueryFormula(:v, :b)
    ]))

    state = SmartsParser("a;b", false)
    and3 = lglowand!(state, parseabc)
    @test and3 == QueryFormula(:and, Set([
        QueryFormula(:v, :a),
        QueryFormula(:v, :b)
    ]))

    state = SmartsParser("", false)
    null = lglowand!(state, parseabc)
    @test null === nothing

    state = SmartsParser("!a!b", false)
    not3 = lglowand!(state, parseabc)
    @test not3 == QueryFormula(:and, Set([
        QueryFormula(:not, QueryFormula(:v, :a)),
        QueryFormula(:not, QueryFormula(:v, :b))
    ]))

    state = SmartsParser("abcdef", false)
    and4 = lghighand!(state, parseabc)
    @test and4 == QueryFormula(:and, Set([
        QueryFormula(:v, :a),
        QueryFormula(:v, :b),
        QueryFormula(:v, :c)
    ]))
    @test state.pos == 4

    state = SmartsParser("a;b&c", false)
    and5 = lglowand!(state, parseabc)
    @test and5 == QueryFormula(:and, Set([
        QueryFormula(:v, :a),
        QueryFormula(:and, Set([
            QueryFormula(:v, :b),
            QueryFormula(:v, :c)
        ]))
    ]))

    state = SmartsParser("a,b&c", false)
    comp1 = lglowand!(state, parseabc)
    @test comp1 == QueryFormula(:or, Set([
        QueryFormula(:v, :a),
        QueryFormula(:and, Set([
            QueryFormula(:v, :b),
            QueryFormula(:v, :c)
        ]))
    ]))

    state = SmartsParser("a,b;c", false)
    comp2 = lglowand!(state, parseabc)
    @test comp2 == QueryFormula(:and, Set([
        QueryFormula(:or, Set([
            QueryFormula(:v, :a),
            QueryFormula(:v, :b)
        ])),
        QueryFormula(:v, :c)
    ]))

    state = SmartsParser("ac;a!b,c;bx", false)
    comp3 = lglowand!(state, parseabc)
    @test comp3 == QueryFormula(:and, Set([
        QueryFormula(:and, Set([
            QueryFormula(:v, :a),
            QueryFormula(:v, :c)
        ])),
        QueryFormula(:or, Set([
            QueryFormula(:and, Set([
                QueryFormula(:v, :a),
                QueryFormula(:not, QueryFormula(:v, :b))
            ])),
            QueryFormula(:v, :c)
        ])),
        QueryFormula(:v, :b)
    ]))
    @test state.pos == 11

end

end # logicalopelator
