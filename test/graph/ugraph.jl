#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.ugraph" begin

@testset "mapugraph" begin
    graph = MapUGraph([1,2,3,4,5], [(1,2), (3,4), (4,5)])
    node = getnode(graph, 4)
    @test typeof(node) <: AbstractNode
    edge = getedge(graph, 3, 4)
    @test typeof(edge) <: AbstractEdge
    @test edge.u == 3
    @test edge.v == 4
    @test getedge(graph, 1, 2) !== Edge(1, 2)
    nbr = neighbors(graph, 4)
    @test length(nbr) == 2
    updatenode!(graph, Node(), 4)
    @test length(neighbors(graph, 4)) == 2
    @test node !== getnode(graph, 4)
    updatenode!(graph, Node(), 6)
    @test length(neighbors(graph, 6)) == 0
    updateedge!(graph, Edge(4, 6), 4)
    @test length(neighbors(graph, 4)) == 3
    @test_throws OperationError updateedge!(graph, Edge(4, 8), 5)
    unlinkedge!(graph, 1, 2)
    @test length(neighbors(graph, 1)) == 0
    @test_throws OperationError unlinkedge!(graph, 7, 8)
    unlinknode!(graph, 4)
    @test length(neighbors(graph, 6)) == 0
    @test_throws OperationError unlinknode!(graph, 8)
    # TODO check pop!(nodekeys)
end

@testset "vectorugraph" begin
    graph = VectorUGraph(6, [(1,2), (3,4), (4,5)])
    node = getnode(graph, 1)
    @test typeof(node) <: AbstractNode
    edge = getedge(graph, 3)
    @test typeof(edge) <: AbstractEdge
    @test edge.u == 4
    @test edge.v == 5
    nbr = neighbors(graph, 6)
    @test length(nbr) == 0
    m = MapUGraph([1,2,3,4,5], [(1,2), (3,4), (4,5)])
    frozen = VectorUGraph{Node,Edge}(m)
    fedge = getedge(frozen, 2, 1)
    @test fedge.u == 1
    @test fedge.v == 2
    @test objectid(m.nodes[5]) == objectid(frozen.nodes[5])
    @test objectid(m.edges[3]) == objectid(frozen.edges[3])
end

end # graph.ugraph
