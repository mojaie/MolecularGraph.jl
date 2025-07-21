
@testset "smarts.logicaloperator" begin

@testset "logicaloperator" begin
    # Test function
    function parseabc!(state, qtree)
        if read(state) == 'a'
            forward!(state)
            return add_qnode!(qtree, qeq(:v, "a"))
        elseif read(state) == 'b'
            forward!(state)
            return add_qnode!(qtree, qeq(:v, "b"))
        elseif read(state) == 'c'
            forward!(state)
            return add_qnode!(qtree, qeq(:v, "c"))
        end
        return 0
    end

    struct QueryTest{T<:Integer,U<:QueryNode} <: QueryTree{T,U}
        graph::SimpleDiGraph{T}
        vprops::Dict{T,U}
    end

    MolecularGraph.lgnot!(state::SMARTSParser, qtree::QueryTest
        ) = lgnot!(state, qtree, parseabc!)

    state = SMARTSParser{Int,QueryTest,QueryTest}("!a")
    qtree = QueryTest(querytree(Tuple{Int,Int}[], QueryNode[])...)
    not1 = lgnot!(state, qtree)
    @test not1 == 2
    @test qtree == QueryTest(querytree([(1, 2)], [qnot(), qeq(:v, "a")])...)

    state = SMARTSParser{Int,QueryAtom,QueryBond}("")
    qtree = QueryTest(querytree(Tuple{Int,Int}[], QueryNode[])...)
    null = lglowand!(state, qtree)
    @test null == 0
    @test qtree == QueryTest(querytree(Tuple{Int,Int}[], QueryNode[])...)

    state = SMARTSParser{Int,QueryAtom,QueryBond}("a&b&cd")
    qtree = QueryTest(querytree(Tuple{Int,Int}[], QueryNode[])...)
    and2 = lghighand!(state, qtree)
    @test and2 == 4
    @test qtree == QueryTest(querytree(
        [(1, 2), (1, 3), (1, 4)],
        [qand(), qeq(:v, "a"), qeq(:v, "b"), qeq(:v, "c")])...)

    state = SMARTSParser{Int,QueryAtom,QueryBond}("!a&b")
    qtree = QueryTest(querytree(Tuple{Int,Int}[], QueryNode[])...)
    not2 = lghighand!(state, qtree)
    @test not2 == 4
    @test qtree == QueryTest(querytree(
        [(1, 2), (2, 3), (1, 4)],
        [qand(), qnot(), qeq(:v, "a"), qeq(:v, "b")])...)

    state = SMARTSParser{Int,QueryAtom,QueryBond}("a,b")
    qtree = QueryTest(querytree(Tuple{Int,Int}[], QueryNode[])...)
    or1 = lgor!(state, qtree)
    @test or1 == 3
    @test qtree == QueryTest(querytree([(1, 2), (1, 3)], [qor(), qeq(:v, "a"), qeq(:v, "b")])...)

    state = SMARTSParser{Int,QueryAtom,QueryBond}("a;b")
    qtree = QueryTest(querytree(Tuple{Int,Int}[], QueryNode[])...)
    and3 = lglowand!(state, qtree)
    @test and3 == 3
    @test qtree == QueryTest(querytree([(1, 2), (1, 3)], [qand(), qeq(:v, "a"), qeq(:v, "b")])...)

    state = SMARTSParser{Int,QueryAtom,QueryBond}("!a!b")
    qtree = QueryTest(querytree(Tuple{Int,Int}[], QueryNode[])...)
    not3 = lglowand!(state, qtree)
    @test not3 == 5
    @test qtree == QueryTest(querytree(
        [(1, 2), (1, 3), (2, 4), (3, 5)],
        [qand(), qnot(), qnot(), qeq(:v, "a"), qeq(:v, "b")])...)

    state = SMARTSParser{Int,QueryAtom,QueryBond}("a&b&c")
    qtree = QueryTest(querytree(Tuple{Int,Int}[], QueryNode[])...)
    and4 = lghighand!(state, qtree)
    @test and4 == 4
    @test qtree == QueryTest(querytree(
        [(1, 2), (1, 3), (1, 4)],
        [qand(), qeq(:v, "a"), qeq(:v, "b"), qeq(:v, "c")])...)
    @test state.pos == 6

    state = SMARTSParser{Int,QueryAtom,QueryBond}("a,b&c")
    qtree = QueryTest(querytree(Tuple{Int,Int}[], QueryNode[])...)
    comp1 = lglowand!(state, qtree)
    @test comp1 == 5
    @test qtree == QueryTest(querytree(
        [(1, 2), (1, 3), (3, 4), (3, 5)],
        [qor(), qeq(:v, "a"), qand(), qeq(:v, "b"), qeq(:v, "c")])...)

    state = SMARTSParser{Int,QueryAtom,QueryBond}("a,b;c")
    qtree = QueryTest(querytree(Tuple{Int,Int}[], QueryNode[])...)
    comp2 = lglowand!(state, qtree)
    @test comp2 == 5
    @test qtree == QueryTest(querytree(
        [(1, 2), (1, 3), (2, 4), (2, 5)],
        [qand(), qor(), qeq(:v, "c"), qeq(:v, "a"), qeq(:v, "b")])...)

    state = SMARTSParser{Int,QueryAtom,QueryBond}("a&c;a!b,c;bx")
    qtree = QueryTest(querytree(Tuple{Int,Int}[], QueryNode[])...)
    comp3 = lglowand!(state, qtree)
    @test comp3 == 11
    @test qtree == QueryTest(querytree(
        [(1, 2), (1, 3), (1, 4), (2, 5), (2, 6), (3, 7), (3, 8), (7, 9), (7, 10), (10, 11)],
        [qand(), qand(), qor(), qeq(:v, "b"), qeq(:v, "a"), qeq(:v, "c"),
        qand(), qeq(:v, "c"), qeq(:v, "a"), qnot(), qeq(:v, "b")])...)
    # (v -> v[1] & v[3] & (v[1] & ~v[2] | v[3]) & v[2]), [(:v, :a), (:v, :b), (:v, :c)]
    @test state.pos == 12
end

end # logicalopelator
