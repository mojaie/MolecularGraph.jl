#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.ugraph" begin

@testset "mapugraph" begin
    graph = MapUDGraph([1,2,3,4,5], [(1,2), (3,4), (4,5)])

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

    nbrnkeys = neighborkeys(graph, 5)
    (nbrn, state) = iterate(nbrnkeys)
    @test nbrn == 4

    nbrnodes = neighbornodes(graph, 5)
    (nbrnode, state) = iterate(nbrnodes)
    @test isa(nbrnode, AbstractNode)

    nbrekeys = neighboredgekeys(graph, 5)
    (nbre, state) = iterate(nbrekeys)
    @test nbre == 3

    nbredges = neighboredges(graph, 5)
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

    newgraph = similarmap(graph)
    @test nodecount(newgraph) == 0
    @test edgecount(newgraph) == 0
    @test typeof(graph) === typeof(newgraph)

    # Nodes and Edges have mutable attribute `attr`
    @test getnode(graph, 1) !== Node()
    newedge = Edge(1, 2)
    @test getedge(graph, 1, 2) !== newedge
    # `connect` function updates connectivity without making new instances
    updated = connect(newedge, 3, 4)
    @test updated.attr === newedge.attr
    @test updated.u == 3
    @test updated.v == 4


    # Conversion to VectorGrpah
    vec = VectorUDGraph{Node,Edge}(graph)
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
    @test node !== getnode(graph, 4)
    # New isolated node
    updatenode!(graph, Node(), 6)
    @test degree(graph, 6) == 0
    # Make new connection
    updateedge!(graph, Edge(4, 6), 4)
    @test degree(graph, 4) == 3
    # Try to connect invalid node
    @test_throws KeyError updateedge!(graph, Edge(4, 8), 5)

    # Delete edge
    unlinkedge!(graph, 1, 2)
    @test degree(graph, 1) == 0
    @test_throws KeyError unlinkedge!(graph, 7, 8)
    # Delete node and adjacent edges
    unlinknode!(graph, 4)
    @test edgecount(graph) == 0
    @test_throws KeyError unlinknode!(graph, 8)
end


@testset "vectorugraph" begin
    graph = VectorUDGraph(6, [(1,2), (3,4), (4,5)])

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

    newgraph = similarmap(graph)
    @test nodecount(newgraph) == 0
    @test edgecount(newgraph) == 0
    @test isa(newgraph, MapUDGraph)
end

end # graph.ugraph
