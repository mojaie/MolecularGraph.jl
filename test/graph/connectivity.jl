#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.connectivity" begin

@testset "connected" begin
    graph0 = mapgraph(1:5, Tuple{Int,Int}[])
    @test issetequal(
        [tuple(c...) for c in connected_components(graph0)],
        [(1,), (2,), (3,), (4,), (5,)]
    )
    graph1 = pathgraph(5)
    @test issetequal(
        [tuple(sort(collect(c))...) for c in connected_components(graph1)],
        [tuple(1:5...)]
    )
    graph2 = mapgraph(1:5, [(1, 2), (2, 3), (4, 5)])
    @test issetequal(
        [tuple(sort(collect(c))...) for c in connected_components(graph2)],
        [(1, 2, 3), (4, 5)]
    )
end


@testset "bridge" begin
    graph1 = pathgraph(5)
    @test issetequal(bridges(graph1), [1, 2, 3, 4])
    graph2 = cyclegraph(5)
    @test isempty(bridges(graph2))
    graph3 = mapgraph(1:8, [
        (1, 2), (2, 3), (1, 3), (3, 4), (4, 5),
        (5, 6), (4, 6), (5, 7), (7, 8), (8, 6)
    ])
    @test issetequal(bridges(graph3), [4])
    graph4 = mapgraph(1:9, [
        (1, 2), (2, 3), (3, 4), (4, 5), (4, 6),
        (3, 7), (7, 8), (8, 9), (9, 1)
    ])
    @test issetequal(bridges(graph4), [3, 4, 5])
    graph5 = mapgraph(1:10, [
        (1, 2), (2, 3), (3, 4), (3, 5), (5, 6),
        (6, 7), (7, 5), (5, 8), (8, 9), (9, 5), (5, 10)
    ])
    @test issetequal(bridges(graph5), [1, 2, 3, 4, 11])
end


@testset "two_edge_connected" begin
    graph1 = pathgraph(5)
    @test isempty(two_edge_connected(graph1))
    graph2 = cyclegraph(5)
    @test issetequal(
        [tuple(sort(collect(c))...) for c in two_edge_connected(graph2)],
        [tuple(1:5...)]
    )
    graph3 = mapgraph(1:8, [
        (1, 2), (2, 3), (1, 3), (3, 4), (4, 5),
        (5, 6), (4, 6), (5, 7), (7, 8), (8, 6)
    ])
    @test issetequal(
        [tuple(sort(collect(c))...) for c in two_edge_connected(graph3)],
        [tuple(1:3...), tuple(4:8...)]
    )
    graph4 = mapgraph(1:9, [
        (1, 2), (2, 3), (3, 4), (4, 5), (4, 6),
        (3, 7), (7, 8), (8, 9), (9, 1)
    ])
    @test issetequal(
        [tuple(sort(collect(c))...) for c in two_edge_connected(graph4)],
        [(1, 2, 3, 7, 8, 9)]
    )
    graph5 = mapgraph(1:10, [
        (1, 2), (2, 3), (3, 4), (3, 5), (5, 6),
        (6, 7), (7, 5), (5, 8), (8, 9), (9, 5), (5, 10)
    ])
    @test issetequal(
        [tuple(sort(collect(c))...) for c in two_edge_connected(graph5)],
        [tuple(5:9...)]
    )
end

end # graph.connectivity
