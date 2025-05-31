#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.multigraph" begin
    p5 = pathgraph(5)
    @test !hasselfloop(p5)
    @test !hasmultiedge(p5)
    @test issimplegraph(p5)
    looped = immutableplaingraph(4, [(1, 2), (2, 3), (3, 1), (3, 3), (3, 4)])
    @test hasselfloop(looped)
    @test !hasmultiedge(looped)
    @test !issimplegraph(looped)
    multi = immutableplaingraph(4, [(1, 2), (2, 3), (3, 1), (3, 4), (3, 4)])
    @test !hasselfloop(multi)
    @test hasmultiedge(multi)
    @test !issimplegraph(multi)
end # graph.multigraph
