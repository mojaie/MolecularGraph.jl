#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.isomorphism.vf2" begin

@testset "subgraph" begin
    g = VectorUDGraph(6, [(1,2), (2,3), (2,4), (2,5), (5,6)])
    h = VectorUDGraph(6, [(1,2), (2,3), (3,4), (3,5), (3,6)])
    @test is_subgraph(h, g)
    @test is_isomorphic(g, h)

    h2 = VectorUDGraph(7, [(1,2), (2,3), (3,4), (3,5), (3,6), (6,7)])
    @test is_subgraph(g, h2)
    @test !is_isomorphic(h2, g)

    h3 = VectorUDGraph(6, [(1,2), (2,3), (3,4), (3,5), (3,6), (6,1)])
    @test !is_subgraph(g, h3)
end

@testset "edgesubgraph" begin
    g = VectorUDGraph(5, [(1,2), (2,3), (2,4), (2,5)])
    h = VectorUDGraph(4, [(1,2), (2,3), (2,4)])
    @test is_edge_subgraph(h, g)

    h2 = VectorUDGraph(5, [(1,2), (2,3), (2,4), (2,5), (5,1)])
    @test is_edge_subgraph(g, h2)

    h3 = VectorUDGraph(5, [(1,2), (2,3), (2,4), (3,5)])
    @test !is_edge_subgraph(g, h3)
end

@testset "deltaY" begin
    tri = cyclegraph(3)
    star = VectorUDGraph(4, [(1,2), (1,3), (1,4)])
    @test is_edge_subgraph(tri, tri)
    @test is_edge_subgraph(star, star)
    @test !is_edge_subgraph(tri, star)
    @test !is_edge_subgraph(star, tri)

    tetra = VectorUDGraph(4, [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)])
    butt = VectorUDGraph(5, [(1,2), (2,3), (1,3), (3,4), (4,5), (5,3)])
    @test !is_edge_subgraph(butt, tetra)

    diam = VectorUDGraph(4, [(1,2), (1,3), (2,3), (2,4), (3,4)])
    eiso = edgeisomorphismiter(tetra, diam, mode=:Subgraph)
    @test length(collect(eiso)) == 24
end

@testset "mandatory" begin
    path = pathgraph(7)
    subp = pathgraph(3)
    eiso = edgeisomorphismiter(path, subp, mode=:Subgraph)
    @test length(collect(eiso)) == 10

    restricted = edgeisomorphismiter(
        path, subp, mode=:Subgraph, mandatory=Dict(3 => 1))
    @test length(collect(restricted)) == 2
end

end # graph.isomorphism.vf2
