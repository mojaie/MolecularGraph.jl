#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "cycle" begin

@testset "circuitrank" begin
    @test circuitrank(SimpleGraph()) == 0
    @test circuitrank(SimpleGraph(3)) == 0
    @test circuitrank(path_graph(5)) == 0
    @test circuitrank(cycle_graph(5)) == 1
    @test circuitrank(complete_graph(5)) == 6
    @test circuitrank(smallgraph(:petersen)) == 6
    twotri = SimpleGraph(Edge.([(1, 2), (2, 3), (3, 1), (4, 5), (5, 6), (6, 4)]))
    @test circuitrank(twotri) == 2
end

@testset "cotree_edges" begin
    @test isempty(cotree_edges(path_graph(5)))
    @test length(cotree_edges(cycle_graph(4))) == 1
    for g in (complete_graph(4), complete_graph(5), smallgraph(:petersen), ladder_graph(5))
        cotree = cotree_edges(g)
        @test length(cotree) == circuitrank(g)
        @test issubset(cotree, Set(edges(g)))
        # the rest is a spanning tree
        tree = SimpleGraph(collect(setdiff(Set(edges(g)), cotree)))
        @test nv(tree) == nv(g)
        @test is_tree(tree)
    end
end

@testset "findmincycle" begin
    c5 = cycle_graph(5)
    cyc = findmincycle(c5, Set([Edge(1, 5)]))
    @test length(cyc) == 5
    @test issetequal(cyc, 1:5)
    # the cycle should contain the edge in S
    g = SimpleGraph(Edge.([(1, 2), (2, 3), (3, 4), (4, 1), (1, 5), (5, 3)]))
    cyc = findmincycle(g, Set([Edge(3, 4)]))
    @test length(cyc) == 4
    @test 3 in cyc && 4 in cyc
end

@testset "mincyclebasis" begin
    p5 = path_graph(5)
    @test isempty(mincyclebasis(p5))
    c5 = cycle_graph(5)
    @test length(mincyclebasis(c5)[1]) == 5
    k44 = complete_bipartite_graph(3, 3)
    @test sum(map(length, mincyclebasis(k44))) == 16
    k5 = complete_graph(5)
    @test sum(map(length, mincyclebasis(k5))) == 18
    cycs = SimpleGraph(Edge.([(1, 2), (2, 3), (3, 1), (3, 4), (4, 5), (5, 6), (6, 4)]))
    @test length(mincyclebasis(cycs)) == 2
    rem_edge!(cycs, 3, 4)
    mcb = mincyclebasis(cycs)
    # connected components are assumed to be in index order
    @test issetequal(mcb[1], [1, 2, 3])
    @test issetequal(mcb[2], [4, 5, 6])
    # https://github.com/mojaie/MolecularGraph.jl/issues/119
    case57 = SimpleGraph(Edge.([
        (1, 2), (1, 15), (1, 16), (1, 17), (2, 3), (3, 4), (3, 15), (4, 5), (4, 6), (4, 18),
        (5, 6), (6, 7), (6, 8), (7, 8), (7, 29), (8, 9), (9, 10), (9, 11), (9, 12), (9, 13),
        (9, 55), (10, 12), (10, 51), (11, 13), (11, 41), (11, 43), (12, 13), (12, 16), (12, 17),
        (13, 14), (13, 15), (13, 49), (14, 15), (14, 46), (15, 45), (18, 19), (19, 20), (20, 21),
        (21, 22), (22, 23), (22, 38), (23, 24), (24, 25), (24, 26), (25, 30), (26, 27), (27, 28),
        (28, 29), (29, 52), (30, 31), (31, 32), (32, 33), (32, 34), (34, 35), (35, 36), (36, 37),
        (36, 40), (37, 38), (37, 39), (38, 44), (38, 48), (38, 49), (39, 57), (40, 56), (41, 42),
        (41, 43), (41, 56), (42, 56), (44, 45), (46, 47), (47, 48), (48, 49), (49, 50), (50, 51),
        (52, 53), (53, 54), (54, 55), (56, 57)]))
    @test sum(length.(mincyclebasis(case57))) == 124
    @test length(mincyclebasis(case57)) == circuitrank(case57)
end

@testset "edgemincyclebasis" begin
    @test isempty(edgemincyclebasis(path_graph(5)))
    ecycs = edgemincyclebasis(cycle_graph(4))
    @test length(ecycs) == 1
    @test issetequal(only(ecycs), edges(cycle_graph(4)))
    k4 = complete_graph(4)
    ecycs = edgemincyclebasis(k4)
    @test length(ecycs) == 3
    @test all(length(c) == 3 for c in ecycs)
    @test all(issubset(c, edges(k4)) for c in ecycs)
    # from vertex cycles
    @test edgemincyclebasis(k4, [[1, 2, 3]]) == [[Edge(1, 2), Edge(2, 3), Edge(1, 3)]]
    @test edgemincyclebasis(k4, Vector{Int}[]) == Vector{Edge{Int}}[]
end

end # cycle
