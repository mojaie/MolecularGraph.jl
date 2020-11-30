#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.inducedsubgraph" begin

@testset "subgraphview" begin
    graph = pathgraph(5)
    subg = nodesubgraph(graph, Set([3, 4, 5]))
    @test_throws KeyError neighbors(subg, 1)[1]
    @test !hasedge(subg, 2, 3)
    @test hasedge(subg, 5, 4)

    # Node and edge key accessor should not be destructive
    nset = nodeset(subg)
    @test issetequal(nset, [3, 4, 5])
    pop!(nset)
    @test issetequal(nodeset(subg), [3, 4, 5])
    eset = edgeset(subg)
    @test issetequal(eset, [3, 4])
    pop!(eset)
    @test issetequal(edgeset(subg), [3, 4])
    @test nodecount(subg) == 3
    @test edgecount(subg) == 2

    graph = completegraph(5)
    nsubg = nodesubgraph(graph, Set([1, 2, 3, 4]))
    @test nodecount(nsubg) == 4
    @test edgecount(nsubg) == 6
    esubg = edgesubgraph(graph, Set([1, 5, 8, 10]))
    @test nodecount(esubg) == 5
    @test edgecount(esubg) == 4
    @test edgesubgraph(graph, [1, 5, 8, 10]) == esubg

    subgsubg = nodesubgraph(nsubg, Set([2, 3, 4]))
    @test issetequal(nodeset(subgsubg), [2, 3, 4])
    @test issetequal(edgeset(subgsubg), [5, 6, 8])
    @test degree(subgsubg, 4) == 2
    newg = plaingraph(subgsubg)
    @test issetequal(nodeset(newg), [1, 2, 3])
    @test issetequal(edgeset(newg), [1, 2, 3])
end

@testset "disubgraphview" begin
    graph = plaindigraph(5, [(1, 2), (2, 3), (3, 4), (4, 5)])
    subg = nodesubgraph(graph, Set([3, 4, 5]))
    @test_throws KeyError neighbors(subg, 1)[1]
    @test !hasedge(subg, 2, 3)
    @test hasedge(subg, 3, 4)
    @test !hasedge(subg, 5, 4)

    @test issetequal(nodeset(subg), [3, 4, 5])
    @test issetequal(edgeset(subg), [3, 4])
    @test nodecount(subg) == 3
    @test edgecount(subg) == 2

    graph = plaindigraph(5,
        [(1, 2), (2, 3), (3, 4), (4, 5), (5, 1),
         (1, 3), (3, 5), (5, 2), (2, 4), (4, 1)])
    nsubg = nodesubgraph(graph, Set([1, 2, 3, 4]))
    @test nodecount(nsubg) == 4
    @test edgecount(nsubg) == 6
    esubg = edgesubgraph(graph, Set([1, 2, 3, 4]))
    @test nodecount(esubg) == 5
    @test edgecount(esubg) == 4

    subgsubg = nodesubgraph(nsubg, Set([2, 3, 4]))
    @test issetequal(nodeset(subgsubg), [2, 3, 4])
    @test issetequal(edgeset(subgsubg), [2, 3, 9])
    @test degree(subgsubg, 4) == 2
    newg = plaindigraph(subgsubg)
    @test issetequal(nodeset(newg), [1, 2, 3])
    @test issetequal(edgeset(newg), [1, 2, 3])
end

end # graph.inducedsubgraph
