#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.component" begin

@testset "connected" begin
    graph0 = MapUGraph(1:5, Tuple{Int,Int}[])
    @test issetequal(
        [tuple(c...) for c in connected_components(graph0)],
        [(1,), (2,), (3,), (4,), (5,)]
    )
    graph1 = MapUGraph(1:5, [(1, 2), (2, 3), (3, 4), (4, 5)])
    @test issetequal(
        [tuple(sort(collect(c))...) for c in connected_components(graph1)],
        [tuple(1:5...)]
    )
    graph2 = MapUGraph(1:5, [(1, 2), (2, 3), (4, 5)])
    @test issetequal(
        [tuple(sort(collect(c))...) for c in connected_components(graph2)],
        [(1, 2, 3), (4, 5)]
    )
end

@testset "two_edge_connected" begin
    graph1 = MapUGraph(1:5, [(1, 2), (2, 3), (3, 4), (4, 5)])
    @test isempty(two_edge_connected(graph1))
    graph2 = MapUGraph(1:5, [(1, 2), (2, 3), (3, 4), (4, 5), (5, 1)])
    @test issetequal(
        [tuple(sort(collect(c))...) for c in two_edge_connected(graph2)],
        [tuple(1:5...)]
    )
    graph3 = MapUGraph(1:8, [
        (1, 2), (2, 3), (1, 3), (3, 4), (4, 5),
        (5, 6), (4, 6), (5, 7), (7, 8), (8, 6)
    ])
    @test issetequal(
        [tuple(sort(collect(c))...) for c in two_edge_connected(graph3)],
        [tuple(1:3...), tuple(4:8...)]
    )
    graph4 = MapUGraph(1:9, [
        (1, 2), (2, 3), (3, 4), (4, 5), (4, 6),
        (3, 7), (7, 8), (8, 9), (9, 1)
    ])
    @test issetequal(
        [tuple(sort(collect(c))...) for c in two_edge_connected(graph4)],
        [(1, 2, 3, 7, 8, 9)]
    )
    graph5 = MapUGraph(1:10, [
        (1, 2), (2, 3), (3, 4), (3, 5), (5, 6),
        (6, 7), (7, 5), (5, 8), (8, 9), (9, 5), (5, 10)
    ])
    @test issetequal(
        [tuple(sort(collect(c))...) for c in two_edge_connected(graph5)],
        [tuple(5:9...)]
    )
end

end # graph.component
