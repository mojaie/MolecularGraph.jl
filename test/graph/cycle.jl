#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.cycle" begin

@testset "mincyclebasis" begin
    p5 = path_graph(5)
    @test isempty(mincyclebasis(p5))
    c5 = cycle_graph(5)
    @test length(mincyclebasis(c5)[1]) == 5
    k44 = complete_bipartite_graph(3, 3)
    @test sum(map(length, mincyclebasis(k44))) == 16
    k5 = complete_graph(5)
    @test sum(map(length, mincyclebasis(k5))) == 18
    cycs = SimpleGraph(Edge.([(1, 2), (2, 3), (3, 1), (3, 4), (4, 5), (5, 6), (6, 4)]))
    @test length(mincyclebasis(cycs)) == 2
end

end # graph.cycle
