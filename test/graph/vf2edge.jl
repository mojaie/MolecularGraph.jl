#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.vf2edge" begin

@testset "edgesubgraph" begin
    g = GVectorUGraph(5, [(1,2), (2,3), (2,4), (2,5)])
    h = GVectorUGraph(4, [(1,2), (2,3), (2,4)])
    @test edge_subgraph_isomorph(g, h)
    h2 = GVectorUGraph(5, [(1,2), (2,3), (2,4), (2,5), (5,1)])
    @test edge_subgraph_isomorph(h2, g)
    h3 = GVectorUGraph(5, [(1,2), (2,3), (2,4), (3,5)])
    @test !edge_subgraph_isomorph(h3, g)
end

@testset "deltaY" begin
    tri = GVectorUGraph(3, [(1,2), (2,3), (1,3)])
    star = GVectorUGraph(4, [(1,2), (1,3), (1,4)])
    @test edge_subgraph_isomorph(tri, tri)
    @test edge_subgraph_isomorph(star, star)
    @test !edge_subgraph_isomorph(tri, star)
    @test !edge_subgraph_isomorph(star, tri)
    tetra = GVectorUGraph(4, [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)])
    butt = GVectorUGraph(5, [(1,2), (2,3), (1,3), (3,4), (4,5), (5,3)])
    @test !edge_subgraph_isomorph(tetra, butt)
    diam = GVectorUGraph(4, [(1,2), (1,3), (2,3), (2,4), (3,4)])
    state = vf2edgesubgraphstate(tetra, diam)
    @test length(collect(isomorphmap!(state))) == 24
end

@testset "mandatory" begin
    path = GVectorUGraph(7, [(1,2), (2,3), (3,4), (4,5), (5,6), (6,7)])
    subp = GVectorUGraph(3, [(1,2), (2,3)])
    state = vf2edgesubgraphstate(path, subp)
    @test length(collect(isomorphmap!(state))) == 10
    state = vf2edgesubgraphstate(path, subp)
    state.mandatory[3] = 1
    @test length(collect(isomorphmap!(state))) == 2
end
end # graph.vf2edge
