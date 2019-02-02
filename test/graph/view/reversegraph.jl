#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.view.reversegraph" begin

@testset "reversegraph" begin
    graph = MapDGraph([1, 2, 3, 4, 5], [(1, 2), (2, 3), (3, 4), (4, 5)])
    rev = ReverseGraph(graph)

    node = getnode(rev, 4)
    @test isa(node, AbstractNode)

    edge = getedge(rev, 3, 4)
    @test isa(edge, DirectedEdge)
    @test edge.source == 4
    @test edge.target == 3
    @test getedge(rev, 3) === edge

    nodes = nodesiter(rev)
    (n, state) = iterate(nodes)
    @test in(n[1], [1, 2, 3, 4, 5])
    @test isa(n[2], AbstractNode)

    edges = edgesiter(rev)
    (e, state) = iterate(edges)
    @test in(e[1], [1, 2, 3, 4])
    @test isa(e[2], DirectedEdge)

    succs = successors(rev, 4)
    (s, state) = iterate(succs)
    @test s[1] == 3
    @test s[2] == 3

    preds = predecessors(rev, 4)
    (p, state) = iterate(preds)
    @test p[1] == 5
    @test p[2] == 4

    @test indegree(rev, 5) == 0
    @test outdegree(rev, 1) == 0

    ntype = nodetype(rev)
    @test ntype <: AbstractNode
    etype = edgetype(rev)
    @test etype <: DirectedEdge

    newgraph = similarmap(rev)
    @test nodecount(newgraph) == 0
    @test edgecount(newgraph) == 0
    @test isa(newgraph, MapDGraph)
end

end # graph.view.reversegraph
