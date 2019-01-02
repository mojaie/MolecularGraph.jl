#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.ugraphview" begin

@testset "ugraphview" begin
    graph = MapUGraph([1, 2, 3, 4, 5], [(1, 2), (2, 3), (3, 4), (4, 5)])
    view = GUGraphView(graph, Set([3, 4, 5]), Set([3, 4]))
    node = getnode(view, 3)
    @test isa(node, Node)
    edge = getedge(view, 3, 4)
    @test isa(edge, Edge)
    @test edge.u == 3
    @test edge.v == 4
    @test degree(view, 3) == 1
    @test_throws OperationError getnode(view, 1)
    @test_throws OperationError getedge(view, 2, 3)
    @test nodecount(view) == 3
    @test edgecount(view) == 2
    nullgraph = similarmap(view)
    @test isa(nullgraph, MapUGraph)
    # TODO check pop!(nodekeys)
end

@testset "nodesubgraph" begin
    graph = MapUGraph([1, 2, 3, 4, 5], [(1, 2), (2, 3), (3, 4), (4, 5)])
    subg = nodesubgraph(graph, [2, 3, 4])
    @test nodecount(subg) == 3
    @test edgecount(subg) == 2
    @test neighborkeys(subg, 2) == [3]
    subgsubg = nodesubgraph(subg, [2, 3])
    @test nodecount(subgsubg) == 2
    @test edgecount(subgsubg) == 1
end

end # graph.ugraphview
