#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.connectivity" begin
    emp = plaingraph(5, Tuple{Int,Int}[])
    @test length(connected_components(emp)) == 5
    @test isempty(bridges(emp))
    @test isempty(cutvertices(emp))
    @test length(biconnected_components(emp)) == 0 # no edges
    @test length(two_edge_connected(emp)) == 5

    pg = pathgraph(5)
    @test length(bridges(pg)) == 4
    @test length(connected_components(pg)) == 1
    @test length(biconnected_components(pg)) == 4
    @test length(two_edge_connected(pg)) == 5

    disconn = plaingraph(5, [(1, 2), (2, 3), (4, 5)])
    @test length(connected_components(disconn)) == 2

    cyc = cyclegraph(5)
    @test isempty(bridges(cyc))
    @test length(biconnected_components(cyc)) == 1
    @test length(two_edge_connected(cyc)) == 1

    fused = plaingraph(8, [
        (1, 2), (2, 3), (1, 3), (3, 4), (4, 5),
        (5, 6), (4, 6), (5, 7), (7, 8), (8, 6)
    ])
    bf = bridges(fused)
    @test length(bf) == 1
    @test bf[1] == 4
    @test length(cutvertices(fused)) == 2
    @test length(biconnected_components(fused)) == 3
    @test length(two_edge_connected(fused)) == 2

    branched = plaingraph(9, [
        (1, 2), (2, 3), (3, 4), (4, 5), (4, 6),
        (3, 7), (7, 8), (8, 9), (9, 1)
    ])
    @test length(bridges(branched)) == 3
    @test length(cutvertices(branched)) == 2

    spiro = plaingraph(10, [
        (1, 2), (2, 3), (3, 4), (3, 5), (5, 6),
        (6, 7), (7, 5), (5, 8), (8, 9), (9, 5), (5, 10)
    ])
    @test length(bridges(spiro)) == 5
    @test length(cutvertices(spiro)) == 3
    @test length(biconnected_components(spiro)) == 7
    @test length(two_edge_connected(spiro)) == 6
end # graph.connectivity
