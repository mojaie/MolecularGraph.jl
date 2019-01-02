#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.vf2edge" begin

@testset "edgesubgraph" begin
    g = VectorUGraph(5, [(1,2), (2,3), (2,4), (2,5)])
    h = VectorUGraph(4, [(1,2), (2,3), (2,4)])
    @test is_edge_subgraph(h, g)
    h2 = VectorUGraph(5, [(1,2), (2,3), (2,4), (2,5), (5,1)])
    @test is_edge_subgraph(g, h2)
    h3 = VectorUGraph(5, [(1,2), (2,3), (2,4), (3,5)])
    @test !is_edge_subgraph(g, h3)
end

@testset "deltaY" begin
    tri = VectorUGraph(3, [(1,2), (2,3), (1,3)])
    star = VectorUGraph(4, [(1,2), (1,3), (1,4)])
    @test is_edge_subgraph(tri, tri)
    @test is_edge_subgraph(star, star)
    @test !is_edge_subgraph(tri, star)
    @test !is_edge_subgraph(star, tri)
    tetra = VectorUGraph(4, [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)])
    butt = VectorUGraph(5, [(1,2), (2,3), (1,3), (3,4), (4,5), (5,3)])
    @test !is_edge_subgraph(butt, tetra)
    diam = VectorUGraph(4, [(1,2), (1,3), (2,3), (2,4), (3,4)])
    @test length(collect(edgeisomorphiter(tetra, diam))) == 24
end

@testset "mandatory" begin
    path = VectorUGraph(7, [(1,2), (2,3), (3,4), (4,5), (5,6), (6,7)])
    subp = VectorUGraph(3, [(1,2), (2,3)])
    @test length(collect(edgeisomorphiter(path, subp))) == 10
    restricted = edgeisomorphiter(
        path, subp, mandatory=Base.ImmutableDict(3 => 1))
    @test length(collect(restricted)) == 2
end

end # graph.vf2edge
