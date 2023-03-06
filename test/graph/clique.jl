#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "clique" begin
    nullg = SimpleGraph()
    @test maximum_clique(nullg) == []

    noedges = SimpleGraph(5)
    @test length(maximum_clique(noedges)) == 1

    g1 = cycle_graph(5)
    @test length(maximum_clique(g1)) == 2

    g2 = SimpleGraph(Edge.([(1, 2), (2, 3), (3, 1), (3, 4), (4, 5)]))
    @test issetequal(maximum_clique(g2), [1, 2, 3])

    g3 = SimpleGraph(
        Edge.([(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4), (4, 5)]))
    @test issetequal(maximum_clique(g3), [1, 2, 3, 4])
end # clique
