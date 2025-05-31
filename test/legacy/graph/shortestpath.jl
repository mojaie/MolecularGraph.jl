#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.shortestpath" begin
    graph = plaingraph(10, [
        (1, 2), (2, 3), (1, 4), (4, 5), (3, 7),
        (7, 4), (7, 8), (8, 9), (9, 10)
    ])
    @test distance(graph, 1, 10) == 5
    @test distance(graph, 1, 1) === nothing
    @test distance(graph, 1, 6) === nothing
    @test length(reachablenodes(graph, 1)) == 8
    @test isreachable(graph, 4, 10)
    @test !isreachable(graph, 6, 1)
    @test eccentricity(graph, 8) == 3
    @test diameter(graph) == 5
    @test count(i->i==Inf, distancematrix(graph)) == 28
    @test shortestpathnodes(graph, 1, 10) == [1, 4, 7, 8, 9, 10]
    @test shortestpathnodes(graph, 5, 8) == [5, 4, 7, 8]
    @test length(longestshortestpathnodes(graph)) == 6

    dgraph = plaindigraph(10, [
        (1, 4), (2, 4), (3, 7), (4, 5), (4, 6),
        (4, 7), (6, 9), (7, 8), (7, 9), (7, 10)
    ])
    @test distance(dgraph, 1, 10) == 3
    @test reversedistance(dgraph, 1, 10) === nothing
    @test length(reachablenodes(dgraph, 4)) == 6
    @test !isreachable(dgraph, 3, 6)
    @test eccentricity(dgraph, 2) == 3
    @test diameter(dgraph) == 3
    @test count(i->i==3.0, distancematrix(dgraph)) == 6
    @test shortestpathnodes(dgraph, 3, 10) == [3, 7, 10]
    @test reverseshortestpathnodes(dgraph, 10, 3) == [10, 7, 3]
end # graph.shortestpath
