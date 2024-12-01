#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "cycle" begin

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
end

end # cycle
