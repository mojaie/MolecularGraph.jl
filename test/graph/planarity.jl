#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

using MolecularGraph.Graph:
    merge!, remove!, planaritytest, outerplanaritytest

@testset "graph.planarity" begin

@testset "merge" begin
    # Empty
    ds1 = Vector{Vector{Int}}[]
    ds2 = Vector{Vector{Int}}[[[1], [2]], [[3], [4]]]
    cotree = Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4)
    @test merge!(ds1, ds2, cotree)
    @test ds1[2][2] == [4]

    # Monochromatic to monochromatic
    cotree = Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4)
    ds1 = Vector{Vector{Int}}[[[1], []]]
    ds2 = Vector{Vector{Int}}[[[2], []], [[3], []], [[4], []]]
    @test merge!(ds1, ds2, cotree)
    @test ds1[2][1] == [2, 3, 4]
    ds1 = Vector{Vector{Int}}[[[1], []], [[2], []], [[3], []]]
    ds2 = Vector{Vector{Int}}[[[4], []]]
    @test merge!(ds1, ds2, cotree)
    @test ds1[4][1] == [4]

    cotree = Dict(11 => 1, 12 => 1, 13 => 1, 14 => 1)
    ds1 = Vector{Vector{Int}}[[[11], []]]
    ds2 = Vector{Vector{Int}}[[[12], []], [[13], []], [[14], []]]
    @test merge!(ds1, ds2, cotree)
    @test length(ds1[1][1]) == 3

    cotree = Dict(1 => 1, 21 => 2, 22 => 2, 3 => 3)
    ds1 = Vector{Vector{Int}}[[[1], []], [[21], []], [[3], []]]
    ds2 = Vector{Vector{Int}}[[[22], []]]
    @test merge!(ds1, ds2, cotree)
    @test ds1[3][1] == [22]
    @test ds1[3][2] == [3]

    cotree = Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4, 5 => 5)
    ds1 = Vector{Vector{Int}}[[[1], []], [[3], []], [[4], []], [[5], []]]
    ds2 = Vector{Vector{Int}}[[[2], []]]
    @test merge!(ds1, ds2, cotree)
    @test ds1[1][1] == [1]
    @test ds1[2][1] == [2]
    @test ds1[2][2] == [3, 4, 5]

    # Merge dichromatic
    cotree = Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4, 5 => 5, 6 => 6)
    ds1 = Vector{Vector{Int}}[[[1], []]]
    ds2 = Vector{Vector{Int}}[[[2, 4], [3]]]
    @test !merge!(ds1, ds2, cotree)
    ds1 = Vector{Vector{Int}}[[[2, 4], [3]]]
    ds2 = Vector{Vector{Int}}[[[5], [6]]]
    @test !merge!(ds1, ds2, cotree)
    ds1 = Vector{Vector{Int}}[[[2, 4], [3]]]
    ds2 = Vector{Vector{Int}}[[[5], []], [[6], []]]
    @test merge!(ds1, ds2, cotree)
    @test ds1[2][1] == [5, 6]

    cotree = Dict(11 => 1, 12 => 1, 21 => 2, 22 => 2, 3 => 3)
    ds1 = Vector{Vector{Int}}[[[11], []], [[12, 21], [3]]]
    ds2 = Vector{Vector{Int}}[[[22], []]]
    @test merge!(ds1, ds2, cotree)
    @test ds1[2][1] == [12, 21, 22]
    @test ds1[2][2] == [3]
    ds1 = Vector{Vector{Int}}[[[12], []], [[11, 3], [21]]]
    ds2 = Vector{Vector{Int}}[[[22], []]]
    @test merge!(ds1, ds2, cotree)
    @test ds1[2][1] == [11, 3]
    @test ds1[2][2] == [21, 22]
end

@testset "remove" begin
    ds = Vector{Vector{Int}}[[[1], []], [[2], [3, 4]]]
    remove!(ds, Set([3]))
    @test ds[2][2] == [4]
    ds = Vector{Vector{Int}}[[[1], []], [[2], [3, 4]]]
    remove!(ds, Set([2, 3]))
    @test ds[2][1] == [4]
    ds = Vector{Vector{Int}}[[[1], []], [[2], [3, 4]]]
    remove!(ds, Set([2, 3, 4]))
    @test length(ds) == 1
end

@testset "isplanar" begin
    @test planaritytest(completebipartite(2,3))
    @test !planaritytest(completebipartite(3,3))
    @test planaritytest(completegraph(4))
    @test !planaritytest(completegraph(5))
    @test planaritytest(laddergraph(20))
    @test planaritytest(circularladder(20))
    @test !planaritytest(moebiusladder(20))
    @test planaritytest(pathgraph(20))
    disconn = plaingraph(11, [
        (1, 2), (2, 3), (3, 4), (6, 7), (6, 8),
        (6, 9), (10, 7), (10, 8), (10, 9), (11, 7), (11, 8), (11, 9)
    ])
    @test !planaritytest(disconn)
end

@testset "isouterplanar" begin
    @test outerplanaritytest(completebipartite(2,2))
    @test !outerplanaritytest(completebipartite(2,3))
    @test outerplanaritytest(completegraph(3))
    @test !outerplanaritytest(completegraph(4))
    @test_broken outerplanaritytest(laddergraph(20))
    @test !outerplanaritytest(circularladder(20))
    @test outerplanaritytest(pathgraph(20))
    disconn = plaingraph(10, [
        (1, 2), (2, 3), (3, 4), (6, 7), (6, 8),
        (6, 9), (10, 7), (10, 8), (10, 9)
    ])
    @test !outerplanaritytest(disconn)
end

end # graph.planarity
