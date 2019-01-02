#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.cycle" begin

@testset "cycleedges" begin
    # TODO: deterministic travarsal
    graph1 = MapUGraph(1:5, [(1, 2), (2, 3), (3, 4), (4, 5)])
    @test isempty(cycleedges(graph1, 1))
    graph2 = MapUGraph(1:5, [(1, 2), (2, 3), (3, 4), (4, 5), (5, 1)])
    @test issetequal(cycleedges(graph2, 1), [5])
    graph3 = MapUGraph(1:8, [
        (1, 2), (2, 3), (1, 3), (3, 4), (4, 5),
        (5, 6), (4, 6), (5, 7), (7, 8), (8, 6)
    ])
    @test issetequal(cycleedges(graph3, 1), [3, 6, 7])
    graph4 = MapUGraph(1:9, [
        (1, 2), (2, 3), (3, 4), (4, 5), (4, 6),
        (3, 7), (7, 8), (8, 9), (9, 1)
    ])
    @test issetequal(cycleedges(graph4, 1), [1])
    graph5 = MapUGraph(1:10, [
        (1, 2), (2, 3), (3, 4), (3, 5), (5, 6),
        (6, 7), (7, 5), (5, 8), (8, 9), (9, 5), (5, 10)
    ])
    @test issetequal(cycleedges(graph5, 1), [5, 8])
end

@testset "minimumcycles" begin
    # TODO: canonical cycle indexing
    graph1 = MapUGraph(1:5, [(1, 2), (2, 3), (3, 4), (4, 5)])
    @test isempty(minimumcycles(graph1))
    graph2 = MapUGraph(1:5, [(1, 2), (2, 3), (3, 4), (4, 5), (5, 1)])
    @test issetequal(minimumcycles(graph2)[1], 1:5)
    graph3 = MapUGraph(1:8, [
        (1, 2), (2, 3), (1, 3), (3, 4), (4, 5),
        (5, 6), (4, 6), (5, 7), (7, 8), (8, 6)
    ])
    graph4 = MapUGraph(1:9, [
        (1, 2), (2, 3), (3, 4), (4, 5), (4, 6),
        (3, 7), (7, 8), (8, 9), (9, 1)
    ])
    graph5 = MapUGraph(1:10, [
        (1, 2), (2, 3), (3, 4), (3, 5), (5, 6),
        (6, 7), (7, 5), (5, 8), (8, 9), (9, 5), (5, 10)
    ])
end

end # graph.cycle
