#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.vf2" begin

@testset "subgraph" begin
    g = GVectorUGraph(6, [(1,2), (2,3), (2,4), (2,5), (5,6)])
    h = GVectorUGraph(6, [(1,2), (2,3), (3,4), (3,5), (3,6)])
    @test subgraph_isomorph(g, h)
    @test is_isomorphic(g, h)
    h2 = GVectorUGraph(7, [(1,2), (2,3), (3,4), (3,5), (3,6), (6,7)])
    @test subgraph_isomorph(h2, g)
    @test !is_isomorphic(h2, g)
    h3 = GVectorUGraph(6, [(1,2), (2,3), (3,4), (3,5), (3,6), (6,1)])
    @test !subgraph_isomorph(h3, g)
end

end # graph.vf2
