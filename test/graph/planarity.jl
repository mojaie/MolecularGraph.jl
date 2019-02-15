#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.planarity" begin

@testset "is_planar" begin
    @test is_planar(bipartitegraph(2,3))
    @test !is_planar(bipartitegraph(3,3))
    @test is_planar(completegraph(4))
    @test !is_planar(completegraph(5))
    @test is_planar(laddergraph(100))
    @test is_planar(circularladder(100))
    @test !is_planar(moebiusladder(10))
end

@testset "is_outerplanar" begin
    @test is_outerplanar(bipartitegraph(2,2))
    @test !is_outerplanar(bipartitegraph(2,3))
    @test is_outerplanar(completegraph(3))
    @test !is_outerplanar(completegraph(4))
    @test is_outerplanar(laddergraph(100))
    @test !is_outerplanar(circularladder(100))
end

end # graph.planarity
