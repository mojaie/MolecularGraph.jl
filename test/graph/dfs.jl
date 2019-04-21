#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.dfs" begin

@testset "circuitrank" begin
    null = plaingraph()
    @test circuitrank(null) == 0
    p5 = pathgraph(5)
    @test circuitrank(p5) == 0
    c5 = cyclegraph(5)
    @test circuitrank(c5) == 1
    lad5 = laddergraph(5)
    @test circuitrank(lad5) == 4
    k5 = completegraph(5)
    @test circuitrank(k5) == 6
    disconn = disjointunion(c5, c5, c5)
    @test circuitrank(disconn) == 3
    # TODO: cached mincycles
    # TODO: subgraphview
end

end # graph.dfs
