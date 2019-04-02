#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.dgraph" begin

@testset "digraph" begin
    graph = digraph(5, [(1,2), (3,4), (4,5)])

    edge = getedge(graph, 2)
    @test edge.source == 3
    @test edge.target == 4
    @test hasedge(graph, 3, 4)
    @test !hasedge(graph, 4, 3)

    succs = outneighbors(graph, 4)
    (s, state) = iterate(succs)
    @test s[1] == 5
    @test s[2] == 3

    preds = outneighbors(graph, 4)
    (p, state) = iterate(preds)
    @test p[1] == 5
    @test p[2] == 3

    nodes = nodesiter(graph)
    (n, state) = iterate(nodes)
    @test isa(n, Tuple)

    edges = edgesiter(graph)
    (e, state) = iterate(edges)
    @test isa(e, Tuple)

    @test indegree(graph, 4) == 1
    @test outdegree(graph, 4) == 1
    @test neighborcount(graph, 4) == 2
    @test degree(graph, 4) == 2

    nkeys = nodekeys(graph)
    @test length(nkeys) == 5
    @test nodecount(graph) == 5

    ekeys = edgekeys(graph)
    @test length(ekeys) == 3
    @test edgecount(graph) == 3

    ntype = nodetype(graph)
    @test ntype <: AbstractNode
    etype = edgetype(graph)
    @test etype <: DirectedEdge

    newgraph = digraph(Node,Arrow)
    @test nodecount(newgraph) == 0
    @test edgecount(newgraph) == 0

    updatenode!(graph, Node())
    @test degree(graph, 6) == 0
    @test_throws DomainError updatenode!(graph, Node(), 8)
    updateedge!(graph, Arrow(4, 6))
    @test outdegree(graph, 4) == 2
    @test_throws DomainError updateedge!(graph, Arrow(4, 8))

    nodesremoved = unlinknodes(graph, (4, 6))
    @test edgecount(nodesremoved) == 1
    edgesremoved = unlinkedges(graph, (1, 2))
    @test nodecount(edgesremoved) == 6
end

end # graph.dgraph
