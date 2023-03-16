
using MolecularGraph:
    forward!, lgnot!, lghighand!, lgor!, lglowand!


@testset "smarts.logicaloperator" begin

@testset "logicaloperator" begin
    # Test function
    function parseabc(state)
        if read(state) == 'a'
            forward!(state)
            return QueryLiteral(:v, :a)
        elseif read(state) == 'b'
            forward!(state)
            return QueryLiteral(:v, :b)
        elseif read(state) == 'c'
            forward!(state)
            return QueryLiteral(:v, :c)
        end
    end

    state = SMARTSParser{SMARTSMolGraph}("!a")
    not1 = lgnot!(state, parseabc)
    @test QueryTruthTable(not1) == QueryTruthTable(v -> ~v[1], [(:v, :a)])

    state = SMARTSParser{SMARTSMolGraph}("a&b")
    and1 = lghighand!(state, parseabc)
    @test QueryTruthTable(and1) == QueryTruthTable(v -> v[1] & v[2], [(:v, :a), (:v, :b)])

    state = SMARTSParser{SMARTSMolGraph}("abc")
    and2 = lghighand!(state, parseabc)
    @test QueryTruthTable(and2) == QueryTruthTable(
        v -> v[1] & v[2] & v[3], [(:v, :a), (:v, :b), (:v, :c)])

    state = SMARTSParser{SMARTSMolGraph}("!ab")
    not2 = lghighand!(state, parseabc)
    @test QueryTruthTable(not2) == QueryTruthTable(v -> ~v[1] & v[2], [(:v, :a), (:v, :b)])

    state = SMARTSParser{SMARTSMolGraph}("a,b")
    or1 = lgor!(state, parseabc)
    @test QueryTruthTable(or1) == QueryTruthTable(v -> v[1] | v[2], [(:v, :a), (:v, :b)])

    state = SMARTSParser{SMARTSMolGraph}("a;b")
    and3 = lglowand!(state, parseabc)
    @test QueryTruthTable(and3) == QueryTruthTable(v -> v[1] & v[2], [(:v, :a), (:v, :b)])

    state = SMARTSParser{SMARTSMolGraph}("")
    null = lglowand!(state, parseabc)
    @test null === nothing

    state = SMARTSParser{SMARTSMolGraph}("!a!b")
    not3 = lglowand!(state, parseabc)
    @test QueryTruthTable(not3) == QueryTruthTable(v -> ~v[1] & ~v[2], [(:v, :a), (:v, :b)])


    state = SMARTSParser{SMARTSMolGraph}("abcdef")
    and4 = lghighand!(state, parseabc)
    @test QueryTruthTable(and4) == QueryTruthTable(
        v -> v[1] & v[2] & v[3], [(:v, :a), (:v, :b), (:v, :c)])
    @test state.pos == 4

    state = SMARTSParser{SMARTSMolGraph}("a,b&c")
    comp1 = lglowand!(state, parseabc)
    @test QueryTruthTable(comp1) == QueryTruthTable(
        v -> v[1] | v[2] & v[3], [(:v, :a), (:v, :b), (:v, :c)])

    state = SMARTSParser{SMARTSMolGraph}("a,b;c")
    comp2 = lglowand!(state, parseabc)
    @test QueryTruthTable(comp2) == QueryTruthTable(
        v -> (v[1] | v[2]) & v[3], [(:v, :a), (:v, :b), (:v, :c)])

    state = SMARTSParser{SMARTSMolGraph}("ac;a!b,c;bx")
    comp3 = lglowand!(state, parseabc)
    @test QueryTruthTable(comp3) == QueryTruthTable(
        v -> v[1] & v[3] & (v[1] & ~v[2] | v[3]) & v[2], [(:v, :a), (:v, :b), (:v, :c)])
    @test state.pos == 11
end

end # logicalopelator
