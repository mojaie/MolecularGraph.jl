#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

using MolecularGraph.MolecularGraphModel: cotree_edges


@testset "graph.cycle" begin

@testset "cycleedges" begin
    # TODO: deterministic travarsal
    graph1 = pathgraph(5)
    @test isempty(cotree_edges(graph1, 1))
    graph2 = cyclegraph(5)
    @test issetequal(cotree_edges(graph2, 1), [5])
    graph3 = vectorgraph(8, [
        (1, 2), (2, 3), (1, 3), (3, 4), (4, 5),
        (5, 6), (4, 6), (5, 7), (7, 8), (8, 6)
    ])
    @test issetequal(cotree_edges(graph3, 1), [3, 6, 7])
    graph4 = vectorgraph(9, [
        (1, 2), (2, 3), (3, 4), (4, 5), (4, 6),
        (3, 7), (7, 8), (8, 9), (9, 1)
    ])
    @test issetequal(cotree_edges(graph4, 1), [1])
    graph5 = vectorgraph(10, [
        (1, 2), (2, 3), (3, 4), (3, 5), (5, 6),
        (6, 7), (7, 5), (5, 8), (8, 9), (9, 5), (5, 10)
    ])
    @test issetequal(cotree_edges(graph5, 1), [5, 8])
end

@testset "minimumcycles" begin
    # TODO: canonical cycle indexing
    # TODO: check that the operation is not destructive
    graph1 = pathgraph(5)
    @test isempty(mincycles(graph1))
    graph2 = cyclegraph(5)
    @test issetequal(mincycles(graph2)[1], 1:5)
    graph3 = vectorgraph(8, [
        (1, 2), (2, 3), (1, 3), (3, 4), (4, 5),
        (5, 6), (4, 6), (5, 7), (7, 8), (8, 6)
    ])
    graph4 = vectorgraph(9, [
        (1, 2), (2, 3), (3, 4), (4, 5), (4, 6),
        (3, 7), (7, 8), (8, 9), (9, 1)
    ])
    graph5 = vectorgraph(10, [
        (1, 2), (2, 3), (3, 4), (3, 5), (5, 6),
        (6, 7), (7, 5), (5, 8), (8, 9), (9, 5), (5, 10)
    ])
end

end # graph.cycle
