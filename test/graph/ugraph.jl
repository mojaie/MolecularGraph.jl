#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.ugraph" begin

@testset "mapgraph" begin
    graph = mapgraph([1,2,3,4,5], [(1,2), (3,4), (4,5)])

    node = getnode(graph, 4)
    @test isa(node, AbstractNode)

    edge = getedge(graph, 3, 4)
    @test isa(edge, UndirectedEdge)
    @test edge.u == 3
    @test edge.v == 4
    @test getedge(graph, 2) === edge

    nodes = nodesiter(graph)
    (n, state) = iterate(nodes)
    @test in(n[1], [1,2,3,4,5])  # The access order is not specified
    @test isa(n[2], AbstractNode)

    edges = edgesiter(graph)
    (e, state) = iterate(edges)
    @test in(e[1], [1,2,3])  # The access order is not specified
    @test isa(e[2], UndirectedEdge)

    nbrs = neighbors(graph, 5)
    (nbr, state) = iterate(nbrs)
    @test nbr[1] == 4
    @test nbr[2] == 3

    nbrnkeys = adjacencies(graph, 5)
    (nbrn, state) = iterate(nbrnkeys)
    @test nbrn == 4

    nbrnodes = adjacentnodes(graph, 5)
    (nbrnode, state) = iterate(nbrnodes)
    @test isa(nbrnode, AbstractNode)

    nbrekeys = incidences(graph, 5)
    (nbre, state) = iterate(nbrekeys)
    @test nbre == 3

    nbredges = incidentedges(graph, 5)
    (nbredge, state) = iterate(nbredges)
    @test isa(nbredge, UndirectedEdge)

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
    @test etype <: UndirectedEdge

    newgraph = mapgraph(Node,Edge)
    @test nodecount(newgraph) == 0
    @test edgecount(newgraph) == 0
    @test typeof(graph) === typeof(newgraph)

    # Nodes and Edges are immutable
    @test getnode(graph, 1) === Node()
    @test getedge(graph, 1, 2) === Edge(1, 2)
    edge = getedge(graph, 1, 2)
    reversed = setnodes(edge, 2, 1)
    @test reversed !== edge
    # Conversion to VectorGrpah
    vec = vectorgraph(graph)
    # Indices are sorted automatically
    vecnodes = nodesiter(vec)
    (n, state) = iterate(vecnodes)
    @test n[1] == 1
    vecedges = edgesiter(vec)
    (e, state) = iterate(vecedges)
    @test e[1] == 1
    # Graph elements are not copied
    @test graph.nodes[5] === vec.nodes[5]
    @test graph.edges[3] === vec.edges[3]


    # Update node
    node = getnode(graph, 4)
    updatenode!(graph, Node(), 4)
    @test degree(graph, 4) == 2
    # New isolated node
    updatenode!(graph, Node(), 6)
    @test degree(graph, 6) == 0
    # Make new connection
    updateedge!(graph, Edge(4, 6), 4)
    @test degree(graph, 4) == 3
    # Try to connect invalid node
    @test_throws KeyError updateedge!(graph, Edge(4, 8), 5)

    """
    # Delete edge
    unlinkedge!(graph, 1, 2)
    @test degree(graph, 1) == 0
    @test_throws KeyError unlinkedge!(graph, 7, 8)
    # Delete node and adjacent edges
    unlinknode!(graph, 4)
    @test edgecount(graph) == 0
    @test_throws KeyError unlinknode!(graph, 8)
    """
end


@testset "vectorugraph" begin
    graph = vectorgraph(6, [(1,2), (3,4), (4,5)])

    node = getnode(graph, 1)
    @test isa(node, AbstractNode)

    edge = getedge(graph, 3)
    @test isa(edge, UndirectedEdge)
    @test edge.u == 4
    @test edge.v == 5

    nodes = nodesiter(graph)
    (n, state) = iterate(nodes)
    @test n[1] == 1
    @test isa(n[2], AbstractNode)

    edges = edgesiter(graph)
    (e, state) = iterate(edges)
    @test e[1] == 1
    @test isa(e[2], UndirectedEdge)

    newgraph = vectorgraph(Node,Edge)
    @test nodecount(newgraph) == 0
    @test edgecount(newgraph) == 0
    @test isa(newgraph, VectorGraph)
end

@testset "cache" begin
    graph = vectorgraph(6, [(1,2), (3,4), (4,5)])
    @cache function hoge(graph)
        return fill("hoge", (nodecount(graph),))
    end
    @test length(hoge(graph)) == 6
    updatenode!(graph, Node())
    @test length(hoge(graph)) == 6
    clearcache!(graph)
    @test length(hoge(graph)) == 7


end

end # graph.ugraph
