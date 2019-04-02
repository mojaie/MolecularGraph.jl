#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.shortestpath" begin

@testset "shortestpath" begin
    graph = vectorgraph(10, [
        (1, 2), (2, 3), (1, 4), (4, 5), (3, 7),
        (7, 4), (7, 8), (8, 9), (9, 10)
    ])
    @test shortestpath(graph, 1, 10) == [1, 4, 7, 8, 9, 10]
    @test shortestpath(graph, 5, 8) == [5, 4, 7, 8]
    @test shortestpath(graph, 1, 1) === nothing
    @test shortestpath(graph, 1, 6) === nothing
end

@testset "dgraph" begin
    graph = digraph(10, [
        (1, 4), (2, 4), (3, 7), (4, 5), (4, 6),
        (4, 7), (6, 9), (7, 8), (7, 9), (7, 10)
    ])
    @test issetequal(reachablenodes(graph, 4), [5, 6, 7, 8, 9, 10])
    plen = pathlength(graph, 1)
    @test plen[6] == 2
    @test plen[10] == 3
end

end # graph.shortestpath
