#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.dag" begin

@testset "dag" begin
    graph = digraph(1:10, [
        (1, 4), (2, 4), (3, 7), (4, 5), (4, 6),
        (4, 7), (6, 9), (7, 8), (7, 9), (7, 10)
    ])
    @test issetequal(ancestors(graph, 7), [1, 2, 3, 4])
    @test issetequal(descendants(graph, 4), [5, 6, 7, 8, 9, 10])
    nodes = topologicalsort(graph)
    @test nodes[4] == 4
    @test nodes[end] in [8, 9, 10]
end

end # graph.shortestpath
