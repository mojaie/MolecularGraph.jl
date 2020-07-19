#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.connectivity" begin
    emp = plaingraph(5, Tuple{Int,Int}[])
    @test length(connectedcomponents(emp)) == 5
    @test isempty(bridges(emp))
    @test isempty(cutvertices(emp))
    @test length(biconnectedcomponents(emp)) == 0 # no edges
    @test length(twoedgeconnectedcomponents(emp)) == 5

    pg = pathgraph(5)
    @test length(bridges(pg)) == 4
    @test length(connectedcomponents(pg)) == 1
    @test length(biconnectedcomponents(pg)) == 4
    @test length(twoedgeconnectedcomponents(pg)) == 5

    disconn = plaingraph(5, [(1, 2), (2, 3), (4, 5)])
    @test length(connectedcomponents(disconn)) == 2

    cyc = cyclegraph(5)
    @test isempty(bridges(cyc))
    @test length(biconnectedcomponents(cyc)) == 1
    @test length(twoedgeconnectedcomponents(cyc)) == 1

    fused = plaingraph(8, [
        (1, 2), (2, 3), (1, 3), (3, 4), (4, 5),
        (5, 6), (4, 6), (5, 7), (7, 8), (8, 6)
    ])
    bf = bridges(fused)
    @test length(bf) == 1
    @test bf[1] == 4
    @test length(cutvertices(fused)) == 2
    @test length(biconnectedcomponents(fused)) == 3
    @test length(twoedgeconnectedcomponents(fused)) == 2

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
    @test length(biconnectedcomponents(spiro)) == 7
    @test length(twoedgeconnectedcomponents(spiro)) == 6
end # graph.connectivity
