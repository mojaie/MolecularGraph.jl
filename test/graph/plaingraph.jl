#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.plaingraph" begin
    g = plaingraph(5, [(1, 2), (3, 4), (4, 5)])
    @test neighbors(g, 4)[3] == 5
    @test getedge(g, 2) == (3, 4)
    @test !hasedge(g, 2, 3)
    @test_throws ErrorException nodeattr(g, 1)
    @test_throws ErrorException edgeattr(g, 1)
    @test_throws ErrorException edgeattr(g, 3, 4)

    @test issetequal(adjacencies(g, 4), [3, 5])
    @test issetequal(incidences(g, 4), [2, 3])
    @test issetequal(nodeset(g), 1:5)
    @test issetequal(edgeset(g), 1:3)
    @test_throws ErrorException nodeattrs(g)
    @test_throws ErrorException edgeattrs(g)

    @test nodecount(g) == 5
    @test edgecount(g) == 3
    @test neighborcount(g, 2) == 1
    @test degree(g, 3) == 1

    addnode!(g)
    @test degree(g, 6) == 0
    addedge!(g, 4, 6)
    @test degree(g, 4) == 3
    @test_throws ErrorException nodeattrtype(g)
    @test_throws ErrorException edgeattrtype(g)

    emp = plaingraph()
    @test nodecount(emp) == 0

    imm = immutableplaingraph(5, [(1, 2), (3, 4), (4, 5)])
    @test_throws MethodError addnode!(imm)
    @test_throws MethodError addedge!(imm, 4, 6)

    # TODO build from attributed graph
    # TODO unlink -> graph.inducedsubgraph
end # graph.plaingraph
