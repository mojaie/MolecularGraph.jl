#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.shortestpath" begin

@testset "shortestpath" begin
    graph = GMapUGraph(1:10, [
        (1, 2), (2, 3), (1, 4), (4, 5), (3, 7),
        (7, 4), (7, 8), (8, 9), (9, 10)
    ])
    @test shortestpath(graph, 1, 10) == [1, 4, 7, 8, 9, 10]
    @test shortestpath(graph, 5, 8) == [5, 4, 7, 8]
    @test shortestpath(graph, 1, 1) === nothing
    @test shortestpath(graph, 1, 6) === nothing
end

end # graph.shortestpath
