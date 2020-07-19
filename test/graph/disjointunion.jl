#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.disjointunion" begin

    G = pathgraph(5)
    H = pathgraph(5)
    I = cyclegraph(5)
    disjointunion!(G, H, I)
    @test nodecount(G) == 15
    @test edgecount(G) == 13
    (g, h, i) = connectedcomponents(G)
    @test length(g) == 5
    @test length(h) == 5
    @test length(i) == 5

    gsub = nodesubgraph(G, Set(1:5))
    U = disjointunion(gsub, H, I)
    @test nodecount(U) == 15
    @test edgecount(U) == 13
    @test length(connectedcomponents(U)) == 3
    @test getunionnode(U, 3, 1) == 11
    @test getunionedge(U, 3, 1) == 9
    @test getsourcenode(U, 15).source == 3
    @test getsourcenode(U, 15).sourcekey == 5
    @test getsourceedge(U, 13).source == 3
    @test getsourceedge(U, 13).sourcekey == 5

    k5 = completegraph(5)
    k5sub = nodesubgraph(k5, Set(1:4))
    U = disjointunion(k5sub, k5sub)
    @test nodecount(U) == 8
    @test edgecount(U) == 12
    @test length(connectedcomponents(U)) == 2

end # graph.disjointunion
