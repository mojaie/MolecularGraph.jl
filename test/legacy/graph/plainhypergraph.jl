#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.plainhypergraph" begin
    g = plainhypergraph(6, [Set([1, 2, 3, 4]), Set([1, 4, 5]), Set([6])])
    @test issetequal(neighbors(g, 4)[2], [1, 4, 5])
    @test issetequal(getedge(g, 1), 1:4)
    @test hasedge(g, 2, 3)
    @test !hasedge(g, 2, 5)
    @test_throws ErrorException nodeattr(g, 1)
    @test_throws ErrorException edgeattr(g, 1)
    @test_throws ErrorException edgeattr(g, 3, 4)

    @test issetequal(nodeset(g), 1:6)
    @test issetequal(edgeset(g), 1:3)
    @test_throws ErrorException nodeattrs(g)
    @test_throws ErrorException edgeattrs(g)

    @test nodecount(g) == 6
    @test edgecount(g) == 3
    @test neighborcount(g, 2) == 1
    @test degree(g, 3) == 1

    addnode!(g)
    @test degree(g, 7) == 0
    addedge!(g, Set([4, 6]))
    @test degree(g, 4) == 3
    @test_throws ErrorException nodeattrtype(g)
    @test_throws ErrorException edgeattrtype(g)

    emp = plainhypergraph()
    @test nodecount(emp) == 0

end # graph.plaingraph
