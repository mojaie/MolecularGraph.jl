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

@testset "isomorphisms" begin
    c6 = cycle_graph(6)
    @test length(collect(isomorphisms(c6, c6))) == 12  # dihedral group
    @test length(collect(isomorphisms(complete_graph(4), complete_graph(4)))) == 24
    @test isempty(isomorphisms(path_graph(4), star_graph(4)))
    @test isempty(isomorphisms(c6, path_graph(6)))
    @test issetequal(
        collect(isomorphisms(path_graph(3), path_graph(3))),
        [Dict(1 => 1, 2 => 2, 3 => 3), Dict(1 => 3, 2 => 2, 3 => 1)])
    for m in isomorphisms(smallgraph(:petersen), smallgraph(:petersen))
        @test length(m) == 10
    end
end

@testset "nodesubgraph_isomorphisms" begin
    c6 = cycle_graph(6)
    p3 = path_graph(3)
    @test length(collect(nodesubgraph_isomorphisms(c6, p3))) == 12
    @test isempty(nodesubgraph_isomorphisms(complete_graph(4), p3))  # P3 is not induced
    for m in nodesubgraph_isomorphisms(c6, p3)
        sub, vmap = induced_subgraph(c6, collect(keys(m)))
        @test ne(sub) == 2
        @test all(has_edge(p3, m[src(e)], m[dst(e)]) for e in edges(c6) if haskey(m, src(e)) && haskey(m, dst(e)))
    end
end

@testset "subgraph_monomorphisms" begin
    p3 = path_graph(3)
    @test length(collect(subgraph_monomorphisms(complete_graph(4), p3))) == 24
    @test length(collect(subgraph_monomorphisms(cycle_graph(6), p3))) == 12
    @test subgraph_is_monomorphic(complete_graph(6), cycle_graph(6))
    @test !subgraph_is_monomorphic(cycle_graph(6), cycle_graph(5))
    @test !subgraph_is_monomorphic(complete_graph(4), cycle_graph(6))
    k4 = complete_graph(4)
    c4 = cycle_graph(4)
    for m in subgraph_monomorphisms(k4, c4)
        rev = Dict(v => k for (k, v) in m)
        @test all(has_edge(k4, rev[src(e)], rev[dst(e)]) for e in edges(c4))
    end
end

@testset "vf2_options" begin
    c6 = cycle_graph(6)
    # vmatch
    label = [1, 2, 1, 2, 1, 2]
    @test length(collect(isomorphisms(c6, c6, vmatch=(a, b) -> label[a] == label[b]))) == 6
    @test is_isomorphic(c6, c6, vmatch=(a, b) -> a == b)
    @test !is_isomorphic(c6, c6, vmatch=(a, b) -> false)
    @test !nodesubgraph_is_isomorphic(c6, path_graph(3), vmatch=(a, b) -> false)
    # ematch
    marked = Edge(1, 2)
    ematch = (e1, e2) -> (e1 == marked) == (e2 == marked)
    @test length(collect(isomorphisms(c6, c6, ematch=ematch))) == 2
    # mandatory and forbidden
    @test length(collect(isomorphisms(c6, c6, mandatory=Dict(1 => 1)))) == 2
    @test all(m[1] == 1 for m in isomorphisms(c6, c6, mandatory=Dict(1 => 1)))
    @test length(collect(isomorphisms(c6, c6, forbidden=Dict(1 => 1)))) == 10
    # timeout
    @test isempty(subgraph_monomorphisms(complete_graph(12), cycle_graph(12), timeout=0))
    iter = subgraph_monomorphisms(complete_graph(12), cycle_graph(12), timeout=0)
    iterate(iter)
    @test iter.timedout
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
    for m in edgesubgraph_isomorphisms(k4, diam)
        @test issubset(keys(m), edges(k4))
        @test issubset(values(m), edges(diam))
        @test length(m) == ne(diam)
    end

    # vertex attribute matching
    p3 = path_graph(3)
    label = [1, 2, 1, 2]
    @test length(collect(edgesubgraph_isomorphisms(path_graph(4), p3))) == 4
    @test length(collect(edgesubgraph_isomorphisms(path_graph(4), p3,
        vmatch=(g, h) -> label[g] == [1, 2, 1][h]))) == 2
    # edge attribute matching
    @test length(collect(edgesubgraph_isomorphisms(path_graph(4), p3,
        ematch=(g, h) -> (g == Edge(1, 2)) == (h == Edge(1, 2))))) == 1
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
