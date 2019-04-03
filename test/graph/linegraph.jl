#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.linegraph" begin
    p5L = linegraph(pathgraph(5))
    @test nodecount(p5L) == 4
    @test edgecount(p5L) == 3
    c5L = linegraph(cyclegraph(5))
    @test nodecount(c5L) == 5
    @test edgecount(c5L) == 5
    lad5L = linegraph(laddergraph(5))
    @test nodecount(lad5L) == 13
    @test edgecount(lad5L) == 22
end # graph.linegraph
