#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "isomorphism_vf2" begin

@testset "subgraph" begin
    null = SimpleGraph()
    @test !is_isomorphic(null, null)
    @test !nodesubgraph_is_isomorphic(null, null)
    @test !edgesubgraph_is_isomorphic(null, null)

    p5 = path_graph(5)
    k5 = complete_graph(5)
    k6 = complete_graph(6)
    @test !nodesubgraph_is_isomorphic(k5, p5)
    @test subgraph_is_monomorphic(k5, p5)
    @test nodesubgraph_is_isomorphic(k6, k5)
    @test !nodesubgraph_is_isomorphic(k5, k6)
    @test is_isomorphic(k5, k5)
    @test !is_isomorphic(k6, k5)

    p7 = path_graph(7)
    disconn1 = SimpleGraph(Edge.([(1,2), (2,3), (4,5), (5,6)]))
    @test nodesubgraph_is_isomorphic(p7, disconn1)
    disconn2 = SimpleGraph(Edge.([(1,2), (3,4), (5,6)]))
    @test !nodesubgraph_is_isomorphic(p7, disconn2)
end

@testset "edgesubgraph" begin

    c4 = cycle_graph(4)
    k4 = complete_graph(4)
    star = star_graph(6)
    @test edgesubgraph_is_isomorphic(k4, c4)
    @test !edgesubgraph_is_isomorphic(k4, star)

    p4 = path_graph(4)
    disconn1 = SimpleGraph(Edge.([(1,2), (3,4)]))
    @test edgesubgraph_is_isomorphic(p4, disconn1)
    disconn2 = SimpleGraph(Edge.([(1,2), (3,4), (4,5)]))
    @test !edgesubgraph_is_isomorphic(p4, disconn2)

    # Delta-Y exchange
    tri = complete_graph(3)
    star = star_graph(4)
    @test edgesubgraph_is_isomorphic(tri, tri)
    @test edgesubgraph_is_isomorphic(star, star)
    @test !edgesubgraph_is_isomorphic(tri, star)
    @test !edgesubgraph_is_isomorphic(star, tri)

    k4 = complete_graph(4)
    diam = smallgraph(:diamond)  # 〼
    butterfly = SimpleGraph(Edge.([(1,2), (2,3), (1,3), (3,4), (4,5), (5,3)]))
    @test !edgesubgraph_is_isomorphic(k4, butterfly)
    matches = edgesubgraph_isomorphisms(k4,diam)
    @test length(collect(matches)) == 24
end

@testset "mandatory" begin
    path = path_graph(7)
    subp = path_graph(3)
    eiso = edgesubgraph_isomorphisms(path, subp)
    @test length(collect(eiso)) == 10

    restricted = edgesubgraph_isomorphisms(path, subp, mandatory=Dict(3 => 1))
    @test length(collect(restricted)) == 2
end

end # isomorphism_vf2
