#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.view.inducedsubgraph" begin

@testset "udsubgraph" begin
    graph = pathgraph(5)
    subg = nodesubgraph(graph, [3, 4, 5])
    @test issetequal(subg.edges, [3, 4])

    node = getnode(subg, 4)
    @test isa(node, AbstractNode)
    @test_throws KeyError getnode(subg, 2)

    edge = getedge(subg, 3, 4)
    @test isa(edge, UndirectedEdge)
    @test edge.u == 3
    @test edge.v == 4
    @test getedge(graph, 3) === edge
    @test_throws KeyError getedge(subg, 1)

    nodes = nodesiter(subg)
    (n, state) = iterate(nodes)
    @test in(n[1], [3, 4, 5])
    @test isa(n[2], AbstractNode)

    edges = edgesiter(subg)
    (e, state) = iterate(edges)
    @test in(e[1], [3, 4])
    @test isa(e[2], UndirectedEdge)

    nbrs = neighbors(subg, 3)
    (nbr, state) = iterate(nbrs)
    @test nbr[1] == 4
    @test nbr[2] == 3
    @test degree(subg, 3) == 1

    nkeys = nodekeys(subg)
    @test nodecount(subg) == 3

    ekeys = edgekeys(subg)
    @test edgecount(subg) == 2

    ntype = nodetype(subg)
    @test ntype <: AbstractNode
    etype = edgetype(subg)
    @test etype <: UndirectedEdge

    newgraph = similarmap(subg)
    @test nodecount(newgraph) == 0
    @test edgecount(newgraph) == 0
    @test isa(newgraph, MapUDGraph)

    # Node and edge key accessor should not be destructive
    pop!(nkeys)
    @test nodecount(subg) == 3
    pop!(ekeys)
    @test edgecount(subg) == 2

    # Edge induced subgraph
    esub = edgesubgraph(graph, [2, 3])
    @test issetequal(esub.nodes, [2, 3, 4])

    # View of view
    subgsubg = nodesubgraph(subg, [4, 5])
    @test nodecount(subgsubg) == 2
    @test edgecount(subgsubg) == 1
    @test_throws KeyError getedge(subgsubg, 3)
    @test degree(subgsubg, 4) == 1
    newsubgsubg = similarmap(subgsubg)
    @test isa(newsubgsubg, MapUDGraph)
end

@testset "dsubgraph" begin
    graph = MapDGraph([1, 2, 3, 4, 5], [(1, 2), (2, 3), (3, 4), (4, 5)])
    subg = nodesubgraph(graph, [3, 4, 5])
    @test issetequal(subg.edges, [3, 4])

    node = getnode(subg, 4)
    @test isa(node, AbstractNode)
    @test_throws KeyError getnode(subg, 2)

    edge = getedge(subg, 3, 4)
    @test isa(edge, DirectedEdge)
    @test edge.source == 3
    @test edge.target == 4
    @test getedge(graph, 3) === edge
    @test_throws KeyError getedge(subg, 1)

    nodes = nodesiter(subg)
    (n, state) = iterate(nodes)
    @test in(n[1], [3, 4, 5])
    @test isa(n[2], AbstractNode)

    edges = edgesiter(subg)
    (e, state) = iterate(edges)
    @test in(e[1], [3, 4])
    @test isa(e[2], DirectedEdge)

    nbrs = neighbors(subg, 3)
    (nbr, state) = iterate(nbrs)
    @test nbr[1] == 4
    @test nbr[2] == 3
    @test degree(subg, 3) == 1

    succs = successors(subg, 3)
    (s, state) = iterate(succs)
    @test s[1] == 4
    @test s[2] == 3
    @test indegree(subg, 3) == 0

    preds = predecessors(subg, 5)
    (p, state) = iterate(preds)
    @test p[1] == 4
    @test p[2] == 4
    @test outdegree(subg, 5) == 0

    ntype = nodetype(subg)
    @test ntype <: AbstractNode
    etype = edgetype(subg)
    @test etype <: DirectedEdge

    newgraph = similarmap(subg)
    @test nodecount(newgraph) == 0
    @test edgecount(newgraph) == 0
    @test isa(newgraph, MapDGraph)

    # Edge induced subgraph
    esub = edgesubgraph(graph, [2, 3])
    @test issetequal(esub.nodes, [2, 3, 4])

    # View of view
    subgsubg = nodesubgraph(subg, [4, 5])
    @test nodecount(subgsubg) == 2
    @test edgecount(subgsubg) == 1
    @test_throws KeyError getedge(subgsubg, 3)
    @test indegree(subgsubg, 4) == 0
    newsubgsubg = similarmap(subgsubg)
    @test isa(newsubgsubg, MapDGraph)
end

end # graph.view.inducedsubgraph
