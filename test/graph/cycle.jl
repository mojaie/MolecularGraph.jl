#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.cycle" begin

@testset "minimumcycles" begin
    # TODO: canonical cycle indexing
    p5 = pathgraph(5)
    @test isempty(edgemincycles(p5))
    c5 = cyclegraph(5)
    @test length(edgemincycles(c5)[1]) == 5
    k44 = completebipartite(3, 3)
    @test sum(map(length, edgemincycles(k44))) == 16
    k5 = completegraph(5)
    @test sum(map(length, edgemincycles(k5))) == 18
end

end # graph.cycle
