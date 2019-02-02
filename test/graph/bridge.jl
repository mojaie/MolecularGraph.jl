#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.bridge" begin

@testset "bridge" begin
    graph1 = pathgraph(5)
    @test issetequal(bridges(graph1, 1), [1, 2, 3, 4])
    graph2 = cyclegraph(5)
    @test isempty(bridges(graph2, 1))
    graph3 = MapUDGraph(1:8, [
        (1, 2), (2, 3), (1, 3), (3, 4), (4, 5),
        (5, 6), (4, 6), (5, 7), (7, 8), (8, 6)
    ])
    @test issetequal(bridges(graph3, 1), [4])
    graph4 = MapUDGraph(1:9, [
        (1, 2), (2, 3), (3, 4), (4, 5), (4, 6),
        (3, 7), (7, 8), (8, 9), (9, 1)
    ])
    @test issetequal(bridges(graph4, 1), [3, 4, 5])
    graph5 = MapUDGraph(1:10, [
        (1, 2), (2, 3), (3, 4), (3, 5), (5, 6),
        (6, 7), (7, 5), (5, 8), (8, 9), (9, 5), (5, 10)
    ])
    @test issetequal(bridges(graph5, 1), [1, 2, 3, 4, 11])
end

end # graph.bridge
