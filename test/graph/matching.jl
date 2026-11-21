#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "matching" begin

@testset "max_matching" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    p7 = path_graph(7)
    @test !is_perfect_matching(p7)
    p8 = path_graph(8)
    @test is_perfect_matching(p8)
    c5 = cycle_graph(5)
    @test !is_perfect_matching(c5)
    c6 = cycle_graph(6)
    @test is_perfect_matching(c6)
    s10 = star_graph(10)
    @test length(max_matching(s10)) == 1
    k7 = complete_graph(7)
    @test length(max_matching(k7)) == 3
    k47 = complete_bipartite_graph(4, 7)
    @test length(max_matching(k47)) == 4
    l8 = ladder_graph(8)
    @test is_perfect_matching(l8)
    dodeca = smallgraph(:dodecahedral)
    @test length(max_matching(dodeca)) == 10
    k20 = complete_graph(20)
    @test length(max_matching(k20)) == 10
    c60 = SimpleGraph(Edge.([
        (1, 2), (1, 5), (1, 9), (2, 3), (2, 12), (3, 4), (3, 15), (4, 5),
        (4, 18), (5, 6), (6, 7), (6, 20), (7, 8), (7, 22), (8, 9), (8, 25),
        (9, 10), (10, 11), (10, 26), (11, 12), (11, 29), (12, 13), (13, 14),
        (13, 30), (14, 15), (14, 33), (15, 16), (16, 17), (16, 34), (17, 18),
        (17, 37), (18, 19), (19, 20), (19, 38), (20, 21), (21, 22), (21, 40),
        (22, 23), (23, 24), (23, 42), (24, 25), (24, 44), (25, 26), (26, 27),
        (27, 28), (27, 45), (28, 29), (28, 47), (29, 30), (30, 31), (31, 32),
        (31, 48), (32, 33), (32, 50), (33, 34), (34, 35), (35, 36), (35, 51),
        (36, 37), (36, 53), (37, 38), (38, 39), (39, 40), (39, 54), (40, 41),
        (41, 42), (41, 55), (42, 43), (43, 44), (43, 57), (44, 45), (45, 46),
        (46, 47), (46, 58), (47, 48), (48, 49), (49, 50), (49, 59), (50, 51),
        (51, 52), (52, 53), (52, 60), (53, 54), (54, 55), (55, 56), (56, 57),
        (56, 60), (57, 58), (58, 59), (59, 60)]))
    @test is_perfect_matching(c60)
    add_vertex!(c60)
    add_edge!(c60, Edge(60, 61))
    @test !is_perfect_matching(c60)
    # global_logger(default_logger)
end

@testset "matching_validity" begin
    is_valid_matching(g, m) = (all(has_edge(g, e) for e in m)
        && allunique(vcat([[src(e), dst(e)] for e in m]...)))
    @test isempty(max_matching(SimpleGraph()))
    @test isempty(max_matching(SimpleGraph(3)))
    for (g, size) in (
            (path_graph(7), 3), (cycle_graph(5), 2), (smallgraph(:petersen), 5),
            (complete_bipartite_graph(4, 7), 4), (smallgraph(:dodecahedral), 10),
            (wheel_graph(8), 4))
        m = max_matching(g)
        @test length(m) == size
        @test is_valid_matching(g, m)
    end
    # blossom (odd cycle) contraction
    g = SimpleGraph(Edge.([(1, 2), (2, 3), (3, 4), (4, 5), (5, 1), (5, 6), (3, 7), (7, 8)]))
    m = max_matching(g)
    @test length(m) == 4
    @test is_valid_matching(g, m)
    # the augmenting path passes through only the root of the contracted blossom
    g = SimpleGraph(Edge.([(1, 3), (1, 5), (1, 6), (2, 3), (2, 4), (3, 4), (3, 7), (4, 7), (6, 7)]))
    m = max_matching(g)
    @test length(m) == 3
    @test is_valid_matching(g, m)
    karate = smallgraph(:karate)
    m = max_matching(karate)
    @test length(m) == 13  # Tutte-Berge bound with S = [1, 2, 3, 4, 33, 34]
    @test is_valid_matching(karate, m)
end

@testset "is_perfect_matching" begin
    p4 = path_graph(4)
    @test is_perfect_matching(p4, Set([Edge(1, 2), Edge(3, 4)]))
    @test !is_perfect_matching(p4, Set([Edge(2, 3)]))
    @test is_perfect_matching(SimpleGraph())
    @test !is_perfect_matching(SimpleGraph(2))
end

end # matching
