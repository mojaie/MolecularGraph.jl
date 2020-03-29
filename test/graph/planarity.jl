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
    @test ds1[4][1] == [11]

    # Monochromatic to dichromatic
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

    # Dichromatic to dichromatic
    cotree = Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4, 5 => 5, 6 => 6)
    ds1 = Vector{Vector{Int}}[[[2, 4], [3]]]
    ds2 = Vector{Vector{Int}}[[[1], []]]
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


@testset "isplanar" begin
    @test planaritytest(completebipartite(2,3))
    @test !planaritytest(completebipartite(3,3))
    @test planaritytest(completegraph(4))
    @test !planaritytest(completegraph(5))
    @test planaritytest(laddergraph(5))
    @test planaritytest(circularladder(5))
    @test !planaritytest(moebiusladder(5))
    # TODO non-biconnected cases (ex. berbell graph)
end

@testset "isouterplanar" begin
    @test outerplanaritytest(completebipartite(2,2))
    @test !outerplanaritytest(completebipartite(2,3))
    @test outerplanaritytest(completegraph(3))
    @test !outerplanaritytest(completegraph(4))
    @test outerplanaritytest(laddergraph(5))
    # @test outerplanaritytest(laddergraph(200))
    @test !outerplanaritytest(circularladder(5))
end

# TODO: test cache

end # graph.planarity
