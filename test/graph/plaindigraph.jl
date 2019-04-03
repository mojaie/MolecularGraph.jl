#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.plaindigraph" begin
    g = plaindigraph(5, [(1, 2), (3, 4), (4, 3), (4, 5)])
    @test neighbors(g, 4)[4] == 5
    @test outneighbors(g, 1)[1] == 2
    @test inneighbors(g, 2)[1] == 1
    @test getedge(g, 2) == (3, 4)
    @test getedge(g, 3) == (4, 3)
    @test !hasedge(g, 2, 1)
    @test_throws ErrorException nodeattr(g, 1)
    @test_throws ErrorException edgeattr(g, 1)
    @test_throws ErrorException edgeattr(g, 3, 4)

    @test issetequal(adjacencies(g, 4), [3, 5])
    @test issetequal(successors(g, 3), [4])
    @test issetequal(predecessors(g, 5), [4])
    @test issetequal(incidences(g, 4), [2, 3, 4])
    @test issetequal(outincidences(g, 4), [3, 4])
    @test issetequal(in_incidences(g, 4), [2])
    @test issetequal(nodeset(g), 1:5)
    @test issetequal(edgeset(g), 1:4)
    @test_throws ErrorException nodeattrs(g)
    @test_throws ErrorException edgeattrs(g)

    @test nodecount(g) == 5
    @test edgecount(g) == 4
    @test neighborcount(g, 2) == 1
    @test degree(g, 4) == 3
    @test outdegree(g, 4) == 2
    @test indegree(g, 4) == 1

    addnode!(g)
    @test degree(g, 6) == 0
    addedge!(g, 4, 6)
    @test outdegree(g, 4) == 3
    @test_throws ErrorException nodeattrtype(g)
    @test_throws ErrorException edgeattrtype(g)

    emp = plaindigraph()
    @test nodecount(emp) == 0

    # TODO build from attributed graph
    # TODO build from subgraph
    # TODO unlink
end # graph.plaindigraph
