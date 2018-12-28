#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.vf2edge" begin

@testset "edgesubgraph" begin
    g = GVectorUGraph(5, [(1,2), (2,3), (2,4), (2,5)])
    h = GVectorUGraph(4, [(1,2), (2,3), (2,4)])
    @test is_edge_subgraph(h, g)
    h2 = GVectorUGraph(5, [(1,2), (2,3), (2,4), (2,5), (5,1)])
    @test is_edge_subgraph(g, h2)
    h3 = GVectorUGraph(5, [(1,2), (2,3), (2,4), (3,5)])
    @test !is_edge_subgraph(g, h3)
end

@testset "deltaY" begin
    tri = GVectorUGraph(3, [(1,2), (2,3), (1,3)])
    star = GVectorUGraph(4, [(1,2), (1,3), (1,4)])
    @test is_edge_subgraph(tri, tri)
    @test is_edge_subgraph(star, star)
    @test !is_edge_subgraph(tri, star)
    @test !is_edge_subgraph(star, tri)
    tetra = GVectorUGraph(4, [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)])
    butt = GVectorUGraph(5, [(1,2), (2,3), (1,3), (3,4), (4,5), (5,3)])
    @test !is_edge_subgraph(butt, tetra)
    diam = GVectorUGraph(4, [(1,2), (1,3), (2,3), (2,4), (3,4)])
    @test length(collect(edgeisomorphiter(tetra, diam))) == 24
end

@testset "mandatory" begin
    path = GVectorUGraph(7, [(1,2), (2,3), (3,4), (4,5), (5,6), (6,7)])
    subp = GVectorUGraph(3, [(1,2), (2,3)])
    @test length(collect(edgeisomorphiter(path, subp))) == 10
    restricted = edgeisomorphiter(
        path, subp, mandatory=Base.ImmutableDict(3 => 1))
    @test length(collect(restricted)) == 2
end

end # graph.vf2edge
