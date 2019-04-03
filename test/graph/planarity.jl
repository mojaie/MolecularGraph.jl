#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

using MolecularGraph.MolecularGraphModel: merge!, remove!

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
    @test isplanar(completebipartite(2,3), fastfilter=false)
    @test !isplanar(completebipartite(3,3), fastfilter=false)
    @test isplanar(completegraph(4), fastfilter=false)
    @test !isplanar(completegraph(5), fastfilter=false)
    @test isplanar(laddergraph(5), fastfilter=false)
    @test isplanar(circularladder(5), fastfilter=false)
    @test !isplanar(moebiusladder(5), fastfilter=false)
    # TODO non-biconnected cases (ex. berbell graph)
end

@testset "isouterplanar" begin
    @test isouterplanar(completebipartite(2,2), fastfilter=false)
    @test !isouterplanar(completebipartite(2,3), fastfilter=false)
    @test isouterplanar(completegraph(3), fastfilter=false)
    @test !isouterplanar(completegraph(4), fastfilter=false)
    @test isouterplanar(laddergraph(5), fastfilter=false)
    # @test isouterplanar(laddergraph(200), fastfilter=false)
    @test !isouterplanar(circularladder(5), fastfilter=false)
end

end # graph.planarity
