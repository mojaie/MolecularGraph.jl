#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "isomorphism_edge" begin

@testset "line_graph_match" begin
    tri = complete_graph(3)
    star = star_graph(4)
    lg, grev, gsh = line_graph(tri)  # grev: 1 => (1, 2), 2 => (1, 3), 3 => (2, 3)
    lh, hrev, hsh = line_graph(star)  # hrev: 1 => (1, 2), 2 => (1, 3), 3 => (1, 4)
    vmatch = lgvmatch(lg, lh, grev, hrev, (a, b) -> true, (a, b) -> true)
    @test all(vmatch(i, j) for i in 1:3, j in 1:3)
    # both end vertices should match (in either direction)
    vmatch = lgvmatch(lg, lh, grev, hrev, (a, b) -> a == b, (a, b) -> true)
    @test [vmatch(i, j) for i in 1:3, j in 1:3] == Bool[1 0 0; 0 1 0; 0 0 0]
    vmatch = lgvmatch(lg, lh, grev, hrev, (a, b) -> true, (a, b) -> a == b)
    @test [vmatch(i, j) for i in 1:3, j in 1:3] == Bool[1 0 0; 0 1 0; 0 0 0]
    # shared vertices of the adjacent edges should match
    ematch = lgematch(lg, lh, gsh, hsh, (a, b) -> a == b)
    @test ematch(Edge(1, 2), Edge(2, 3))  # tri vertex 1, star vertex 1
    @test !ematch(Edge(2, 3), Edge(1, 2))  # tri vertex 3, star vertex 1

    # vertex attribute vectors of the line graph (order independent)
    vvec = lgvmatchvecgen(grev, v -> v, (u, v) -> 0)
    @test [vvec(i) for i in 1:3] == [1 * 120 + 2, 1 * 120 + 3, 2 * 120 + 3]
    vvec = lgvmatchvecgen(grev, v -> 4 - v, (u, v) -> 1)
    @test [vvec(i) for i in 1:3] == [14400 + 2 * 120 + 3, 14400 + 1 * 120 + 3, 14400 + 1 * 120 + 2]
    evec = lgematchvecgen(lg, gsh, v -> 10v)
    @test [evec(i, j) for i in 1:3, j in 1:3] == [0 10 20; 10 0 30; 20 30 0]
    evec = lgematchvecgen(lh, hsh, v -> 10v)
    @test evec(1, 3) == 10
end

@testset "delta_y_edges" begin
    @test isempty(delta_edges(SimpleGraph()))
    @test isempty(delta_edges(cycle_graph(4)))
    tris = delta_edges(complete_graph(3))
    @test !isempty(tris)
    @test all(issetequal(t, [Edge(1, 2), Edge(1, 3), Edge(2, 3)]) for t in tris)
    @test length(unique(Set.(delta_edges(complete_graph(4))))) == 4
    # non-uniform vertex attributes
    @test isempty(delta_edges(complete_graph(3), (a, b) -> a + b != 3))

    @test isempty(y_edges(SimpleGraph()))
    @test isempty(y_edges(path_graph(3)))
    ys = y_edges(star_graph(4))
    @test length(ys) == 1
    @test issetequal(only(ys), [Edge(1, 2), Edge(1, 3), Edge(1, 4)])
    @test length(y_edges(star_graph(5))) == 4  # C(4, 3)
    @test isempty(y_edges(star_graph(4), (a, b) -> b != 4))
end

@testset "delta_y_test" begin
    tri = complete_graph(3)
    star = star_graph(4)
    gd, gy, hd, hy = delta_edges(tri), y_edges(tri), delta_edges(star), y_edges(star)
    # Δ-Y mismatch
    mapping = Dict(Edge(1, 2) => Edge(1, 2), Edge(1, 3) => Edge(1, 3), Edge(2, 3) => Edge(1, 4))
    @test !delta_y_test(mapping, gd, gy, hd, hy)
    revmap = Dict(v => k for (k, v) in mapping)
    @test !delta_y_test(revmap, hd, hy, gd, gy)
    # two edges can be matched
    partial = Dict(Edge(1, 2) => Edge(1, 2), Edge(1, 3) => Edge(1, 3))
    @test delta_y_test(partial, gd, gy, hd, hy)
    @test delta_y_test(Dict(v => k for (k, v) in partial), hd, hy, gd, gy)

    # one of the triangle edges is removed
    corrected = copy(mapping)
    delta_y_correction!(corrected, gd, gy, hd, hy)
    @test length(corrected) == 2
    @test delta_y_test(corrected, gd, gy, hd, hy)
    corrected = copy(revmap)
    delta_y_correction!(corrected, hd, hy, gd, gy)
    @test length(corrected) == 2
    @test delta_y_test(corrected, hd, hy, gd, gy)
    # no change
    corrected = copy(partial)
    delta_y_correction!(corrected, gd, gy, hd, hy)
    @test corrected == partial
end

end # isomorphism_edge
