#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.dgraph" begin

@testset "mapdgraph" begin
    graph = MapDGraph([1,2,3,4,5], [(1,2), (3,4), (4,5)])
    node = getnode(graph, 4)
    @test typeof(node) <: AbstractNode
    edge = getedge(graph, 3, 4)
    @test typeof(edge) <: AbstractDirectedEdge
    @test edge.source == 3
    @test edge.target == 4
    @test getedge(graph, 1, 2) !== Arrow(1, 2)
    succs = successors(graph, 4)
    @test length(succs) == 1
    preds = predecessors(graph, 4)
    @test length(preds) == 1
    updatenode!(graph, Node(), 4)
    @test length(successors(graph, 4)) == 1
    @test node !== getnode(graph, 4)
    updatenode!(graph, Node(), 6)
    @test length(successors(graph, 6)) == 0
    @test length(predecessors(graph, 6)) == 0
    updateedge!(graph, Arrow(4, 6), 4)
    @test length(successors(graph, 4)) == 2
    @test_throws OperationError updateedge!(graph, Arrow(4, 8), 5)
    unlinkedge!(graph, 1, 2)
    @test length(successors(graph, 1)) == 0
    @test_throws OperationError unlinkedge!(graph, 7, 8)
    unlinknode!(graph, 4)
    @test length(successors(graph, 6)) == 0
    @test_throws OperationError unlinknode!(graph, 8)
    # TODO check pop!(nodekeys)
end

end # graph.dgraph
