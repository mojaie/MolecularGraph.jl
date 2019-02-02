#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.dgraph" begin

@testset "mapdgraph" begin
    graph = MapDGraph([1,2,3,4,5], [(1,2), (3,4), (4,5)])

    node = getnode(graph, 4)
    @test isa(node, AbstractNode)

    edge = getedge(graph, 3, 4)
    @test isa(edge, DirectedEdge)
    @test edge.source == 3
    @test edge.target == 4
    @test getedge(graph, 2) === edge

    nodes = nodesiter(graph)
    (n, state) = iterate(nodes)
    @test in(n[1], [1,2,3,4,5])  # The access order is not specified
    @test isa(n[2], AbstractNode)

    edges = edgesiter(graph)
    (e, state) = iterate(edges)
    @test in(e[1], [1,2,3])  # The access order is not specified
    @test isa(e[2], DirectedEdge)

    succs = successors(graph, 4)
    (s, state) = iterate(succs)
    @test s[1] == 5
    @test s[2] == 3

    preds = successors(graph, 4)
    (p, state) = iterate(preds)
    @test p[1] == 5
    @test p[2] == 3

    nbrs = neighbors(graph, 5)
    (nbr, state) = iterate(nbrs)
    @test nbr[1] == 4
    @test nbr[2] == 3

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

    newgraph = similarmap(graph)
    @test nodecount(newgraph) == 0
    @test edgecount(newgraph) == 0
    @test typeof(graph) === typeof(newgraph)

    # Update node
    node = getnode(graph, 4)
    updatenode!(graph, Node(), 4)
    @test degree(graph, 4) == 2
    @test node !== getnode(graph, 4)
    # New isolated node
    updatenode!(graph, Node(), 6)
    @test degree(graph, 6) == 0
    # Make new connection
    updateedge!(graph, Arrow(4, 6), 4)
    @test outdegree(graph, 4) == 2
    # Try to connect invalid node
    @test_throws KeyError updateedge!(graph, Arrow(4, 8), 5)

    # Delete edge
    unlinkedge!(graph, 1, 2)
    @test degree(graph, 1) == 0
    @test_throws KeyError unlinkedge!(graph, 7, 8)
    # Delete node and adjacent edges
    unlinknode!(graph, 4)
    @test edgecount(graph) == 0
    @test_throws KeyError unlinknode!(graph, 8)
end

end # graph.dgraph
