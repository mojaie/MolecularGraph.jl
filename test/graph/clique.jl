#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.clique" begin
    nullg = plaingraph()
    @test issetequal(maximumclique(nullg), [])

    noedges = plaingraph(5, Tuple{Int,Int}[])
    @test length(maximumclique(noedges)) == 1

    g1 = cyclegraph(5)
    @test length(maximumclique(g1)) == 2

    g2 = plaingraph(5, [(1, 2), (2, 3), (3, 1), (3, 4), (4, 5)])
    @test issetequal(maximumclique(g2), [1, 2, 3])

    g3 = plaingraph(
        5, [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4), (4, 5)])
    @test issetequal(maximumclique(g3), [1, 2, 3, 4])
end # graph.clique
