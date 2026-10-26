#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "operator" begin

@testset "edge_rank" begin
    @test isempty(edge_rank(SimpleGraph(3)))
    c4 = cycle_graph(4)
    ernk = edge_rank(c4)
    @test ernk == Dict(e => i for (i, e) in enumerate(edges(c4)))
    @test edge_rank(ernk, 1, 2) == 1
    @test edge_rank(ernk, 4, 1) == edge_rank(ernk, 1, 4) == 2
    @test edge_rank(ernk, 3, 4) == 4
    @test_throws KeyError edge_rank(ernk, 1, 3)
end

@testset "induced_subgraph_edges" begin
    k4 = complete_graph(4)
    @test issetequal(induced_subgraph_edges(k4, [1, 3, 4]), [Edge(1, 3), Edge(1, 4), Edge(3, 4)])
    @test issetequal(induced_subgraph_edges(path_graph(5), [4, 2, 3]), [Edge(2, 3), Edge(3, 4)])
    @test isempty(induced_subgraph_edges(path_graph(5), [1, 3, 5]))
    @test isempty(induced_subgraph_edges(path_graph(5), Int[]))
end

@testset "disjoint_union" begin
    g = path_graph(3)
    disjoint_union!(g, cycle_graph(3))
    @test nv(g) == 6
    @test ne(g) == 5
    @test issetequal(edges(g), Edge.([(1, 2), (2, 3), (4, 5), (4, 6), (5, 6)]))
    g = path_graph(3)
    disjoint_union!(g, SimpleGraph(2))
    @test nv(g) == 5
    @test ne(g) == 2
    g = SimpleGraph()
    disjoint_union!(g, path_graph(2))
    @test nv(g) == 2
    @test ne(g) == 1

    g = path_graph(2)
    u = disjoint_union(g, complete_graph(3))
    @test nv(g) == 2  # not changed
    @test nv(u) == 5
    @test ne(u) == 4
    @test length(connected_components(u)) == 2
    d = disjoint_union(SimpleDiGraph(Edge.([(1, 2)])), SimpleDiGraph(Edge.([(2, 1)])))
    @test nv(d) == 4
    @test issetequal(edges(d), Edge.([(1, 2), (4, 3)]))
end

@testset "modularproduct" begin
    nul = SimpleGraph{Int}()
    prod, isconn = modular_product(nul, nul)
    @test nv(prod) == 0
    @test ne(prod) == 0
    @test isempty(isconn)
    @test eltype(prod) === Int
    g = path_graph(3)
    prod, isconn = modular_product(g, g)
    @test nv(prod) == 9
    @test ne(prod) == 10
    @test count(values(isconn)) == 8
    # vertex attribute matching
    prod, isconn = modular_product(g, g, vmatch=(g1, h1) -> g1 == h1)
    @test ne(prod) == 3
    @test isconn == Dict(Edge(1, 5) => true, Edge(5, 9) => true, Edge(1, 9) => false)
    # custom edge filter
    prod, isconn = modular_product(g, g, edgefilter=(g1, g2, h1, h2) -> false)
    @test nv(prod) == 9
    @test ne(prod) == 0
    # different graphs
    # (P3 has 2 adjacent and 1 non-adjacent node pairs, K3 has 3 adjacent pairs)
    @test ne(modular_product(path_graph(3), complete_graph(3))[1]) == 12
end

@testset "line_graph" begin
    nul, revmap, shared = line_graph(SimpleGraph{Int}())
    @test nv(nul) == 0
    @test ne(nul) == 0
    @test length(revmap) == 0
    @test length(shared) == 0
    nul, revmap, shared = line_graph(SimpleGraph(10))
    @test nv(nul) == 0
    @test ne(nul) == 0
    p5L, revmap, shared = line_graph(path_graph(5))
    @test nv(p5L) == 4
    @test ne(p5L) == 3
    @test length(revmap) == 4
    @test length(shared) == 3
    c5L, revmap, shared = line_graph(cycle_graph(5))
    @test nv(c5L) == 5
    @test ne(c5L) == 5
    lad5L, revmap, shared = line_graph(ladder_graph(5))
    @test nv(lad5L) == 13
    @test ne(lad5L) == 22
    k3L, revmap, shared = line_graph(complete_graph(3))
    @test revmap == Dict(1 => Edge(1, 2), 2 => Edge(1, 3), 3 => Edge(2, 3))
    @test shared == Dict(Edge(1, 2) => 1, Edge(1, 3) => 2, Edge(2, 3) => 3)
    s4L, revmap, shared = line_graph(star_graph(4))
    @test ne(s4L) == 3
    @test all(==(1), values(shared))
end

end
