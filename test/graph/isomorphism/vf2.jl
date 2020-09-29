#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.isomorphism.vf2" begin

@testset "subgraph" begin
    null = immutableplaingraph()
    @test !isexactmatch(null, null)
    @test !issubgraphmatch(null, null)

    p5 = pathgraph(5)
    k5 = completegraph(5)
    k6 = completegraph(6)
    @test !issubgraphmatch(p5, k5)
    @test issubgraphmatch(k6, k5)
    @test !issubgraphmatch(k5, k6)
    @test isexactmatch(k5, k5)
    @test !isexactmatch(k6, k5)

    p7 = pathgraph(7)
    disconn1 = immutableplaingraph(6, [(1,2), (2,3), (4,5), (5,6)])
    @test issubgraphmatch(p7, disconn1)
    disconn2 = immutableplaingraph(6, [(1,2), (3,4), (5,6)])
    @test !issubgraphmatch(p7, disconn2)

    cd6 = circularladder(6)
    mb6 = moebiusladder(6)
    @test !isexactmatch(cd6, mb6)
    # TODO: relabeled mb6
    # @test isexactmatch(rmb6, mb6)
end

@testset "edgesubgraph" begin
    null = immutableplaingraph()
    @test !isedgesubgraphmatch(null, null)

    c4 = cyclegraph(4)
    k4 = completegraph(4)
    star = immutableplaingraph(6, [(1,2), (1,3), (1,4), (1,5), (1,6)])
    @test isedgesubgraphmatch(k4, c4)
    @test !isedgesubgraphmatch(k4, star)

    p4 = pathgraph(4)
    disconn1 = immutableplaingraph(4, [(1,2), (3,4)])
    @test isedgesubgraphmatch(p4, disconn1)
    disconn2 = immutableplaingraph(5, [(1,2), (3,4), (4,5)])
    @test !isedgesubgraphmatch(p4, disconn2)

    # Delta-Y exchange
    tri = cyclegraph(3)
    star = immutableplaingraph(4, [(1,2), (1,3), (1,4)])
    @test isedgesubgraphmatch(tri, tri)
    @test isedgesubgraphmatch(star, star)
    @test !isedgesubgraphmatch(tri, star)
    @test !isedgesubgraphmatch(star, tri)

    k4 = completegraph(4)
    〼 = immutableplaingraph(4, [(1,2), (2,3), (1,3), (3,4), (4,2)])
    butterfly = immutableplaingraph(
        5, [(1,2), (2,3), (1,3), (3,4), (4,5), (5,3)])
    @test !isedgesubgraphmatch(k4, butterfly)
    matches = edgesubgraphmatches(k4, 〼)
    @test length(collect(matches)) == 24
end

@testset "mandatory" begin
    path = pathgraph(7)
    subp = pathgraph(3)
    eiso = edgesubgraphmatches(path, subp)
    @test length(collect(eiso)) == 10

    restricted = edgesubgraphmatches(path, subp, mandatory=Dict(3 => 1))
    @test length(collect(restricted)) == 2
end

end # graph.isomorphism.vf2
