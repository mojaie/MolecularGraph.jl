#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.isomorphism.vf2" begin

@testset "subgraph" begin
    null = immutableplaingraph()
    @test !is_isomorphic(null, null)
    @test !nodesubgraph_is_isomorphic(null, null)

    p5 = pathgraph(5)
    k5 = completegraph(5)
    k6 = completegraph(6)
    @test !nodesubgraph_is_isomorphic(k5, p5)
    @test subgraph_is_monomorphic(k5, p5)
    @test nodesubgraph_is_isomorphic(k6, k5)
    @test !nodesubgraph_is_isomorphic(k5, k6)
    @test is_isomorphic(k5, k5)
    @test !is_isomorphic(k6, k5)

    p7 = pathgraph(7)
    disconn1 = immutableplaingraph(6, [(1,2), (2,3), (4,5), (5,6)])
    @test nodesubgraph_is_isomorphic(p7, disconn1)
    disconn2 = immutableplaingraph(6, [(1,2), (3,4), (5,6)])
    @test !nodesubgraph_is_isomorphic(p7, disconn2)

    cd6 = circularladder(6)
    mb6 = moebiusladder(6)
    @test !is_isomorphic(cd6, mb6)
    # TODO: relabeled mb6
    # @test isexactmatch(rmb6, mb6)
end

@testset "edgesubgraph" begin
    null = immutableplaingraph()
    @test !edgesubgraph_is_isomorphic(null, null)

    c4 = cyclegraph(4)
    k4 = completegraph(4)
    star = immutableplaingraph(6, [(1,2), (1,3), (1,4), (1,5), (1,6)])
    @test edgesubgraph_is_isomorphic(k4, c4)
    @test !edgesubgraph_is_isomorphic(k4, star)

    p4 = pathgraph(4)
    disconn1 = immutableplaingraph(4, [(1,2), (3,4)])
    @test edgesubgraph_is_isomorphic(p4, disconn1)
    disconn2 = immutableplaingraph(5, [(1,2), (3,4), (4,5)])
    @test !edgesubgraph_is_isomorphic(p4, disconn2)

    # Delta-Y exchange
    tri = cyclegraph(3)
    star = immutableplaingraph(4, [(1,2), (1,3), (1,4)])
    @test edgesubgraph_is_isomorphic(tri, tri)
    @test edgesubgraph_is_isomorphic(star, star)
    @test !edgesubgraph_is_isomorphic(tri, star)
    @test !edgesubgraph_is_isomorphic(star, tri)

    k4 = completegraph(4)
    〼 = immutableplaingraph(4, [(1,2), (2,3), (1,3), (3,4), (4,2)])
    butterfly = immutableplaingraph(
        5, [(1,2), (2,3), (1,3), (3,4), (4,5), (5,3)])
    @test !edgesubgraph_is_isomorphic(k4, butterfly)
    matches = edgesubgraph_isomorphisms(k4, 〼)
    @test length(collect(matches)) == 24
end

@testset "mandatory" begin
    path = pathgraph(7)
    subp = pathgraph(3)
    eiso = edgesubgraph_isomorphisms(path, subp)
    @test length(collect(eiso)) == 10

    restricted = edgesubgraph_isomorphisms(path, subp, mandatory=Dict(3 => 1))
    @test length(collect(restricted)) == 2
end

end # graph.isomorphism.vf2
