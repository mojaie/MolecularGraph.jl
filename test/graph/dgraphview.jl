#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.dgraphview" begin

@testset "dgraphview" begin
    graph = GMapDGraph([1, 2, 3, 4, 5], [(1, 2), (2, 3), (3, 4), (4, 5)])
    view = GDGraphView(graph, Set([3, 4, 5]), Set([3, 4]))
    node = getnode(view, 3)
    @test isa(node, Node)
    edge = getedge(view, 3, 4)
    @test isa(edge, Arrow)
    @test edge.source == 3
    @test edge.target == 4
    @test indegree(view, 4) == 1
    @test outdegree(view, 4) == 1
    @test_throws OperationError getnode(view, 1)
    @test_throws OperationError getedge(view, 2, 3)
    @test nodecount(view) == 3
    @test edgecount(view) == 2
    nullgraph = similarmap(view)
    @test isa(nullgraph, MapDGraph)
    # TODO check pop!(nodekeys)
end

@testset "nodesubgraph" begin
    graph = GMapDGraph([1, 2, 3, 4, 5], [(1, 2), (2, 3), (3, 4), (4, 5)])
    subg = nodesubgraph(graph, [2, 3, 4])
    @test nodecount(subg) == 3
    @test edgecount(subg) == 2
    @test succkeys(subg, 2) == [3]
    @test predkeys(subg, 4) == [3]
    subgsubg = nodesubgraph(subg, [2, 3])
    @test nodecount(subgsubg) == 2
    @test edgecount(subgsubg) == 1
end

end # graph.ugraphview
