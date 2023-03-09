#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

using MolecularGraph: merge_ds!

@testset "graph.planarity" begin

@testset "merge" begin
    cotree = Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4)
    # Any branches can be merged to empty trunk
    ds1 = Vector{Vector{Int}}[]
    ds2 = Vector{Vector{Int}}[[[1], []], [[2, 4], [3]]]
    @test merge_ds!(ds1, ds2, cotree)
    @test ds1 == [[[1], []], [[2, 4], [3]]]

    # ds1 edge has a smaller low(e) than ds2 edges, so ds2 edges are fused
    ds1 = Vector{Vector{Int}}[[[1], []]]
    ds2 = Vector{Vector{Int}}[[[2], []], [[3], []], [[4], []]]
    @test merge_ds!(ds1, ds2, cotree)
    @test ds1 == [[[1], []], [[2, 3, 4], []]]

    # ds1 edges have smaller low(e) than ds2 edge, so ds1 edges are not fused
    ds1 = Vector{Vector{Int}}[[[1], []], [[2], []], [[3], []]]
    ds2 = Vector{Vector{Int}}[[[4], []]]
    @test merge_ds!(ds1, ds2, cotree)
    @test ds1 == [[[1], []], [[2], []], [[3], []], [[4], []]]

    # ds1 edges with greater low(e) than ds2 are fused. d2 is set to the opposite side
    ds1 = Vector{Vector{Int}}[[[1], []], [[3], []], [[4], []]]
    ds2 = Vector{Vector{Int}}[[[2], []]]
    @test merge_ds!(ds1, ds2, cotree)
    @test ds1 == [[[1], []], [[3, 4], [2]]]

    # ds2 edges are fused and set to the opposite side of the ds1 cell that have greater low(e)
    ds1 = Vector{Vector{Int}}[[[1], []], [[3], []]]
    ds2 = Vector{Vector{Int}}[[[2], []], [[4], []]]
    @test merge_ds!(ds1, ds2, cotree)
    @test ds1 == [[[1], []], [[3], [2, 4]]]

    # dichromatic branches can not be merged
    ds1 = Vector{Vector{Int}}[[[1], []]]
    ds2 = Vector{Vector{Int}}[[[2], []], [[4], [3]]]
    @test !merge_ds!(ds1, ds2, cotree)

    # ds2 edge with smaller low(e) can not be merged to the dichromatic ds1 cell
    ds1 = Vector{Vector{Int}}[[[1], []], [[4], [3]]]
    ds2 = Vector{Vector{Int}}[[[2], []]]
    @test !merge_ds!(ds1, ds2, cotree)

    # if low(e) of ds2 bottom is the same as that of ds1 bottom, they are merged into the same cell
    cotree = Dict(11 => 1, 12 => 1, 2 => 2, 3 => 3)
    ds1 = Vector{Vector{Int}}[[[11], []]]
    ds2 = Vector{Vector{Int}}[[[12], []], [[2], []], [[3], []]]
    @test merge_ds!(ds1, ds2, cotree)
    @test ds1 == [[[11, 12], []], [[2, 3], []]]

    # only the top cell in ds1 has greater low(e) than ds2. d2 is set to the opposite side
    cotree = Dict(1 => 1, 21 => 2, 22 => 2, 3 => 3)
    ds1 = Vector{Vector{Int}}[[[1], []], [[21], []], [[3], []]]
    ds2 = Vector{Vector{Int}}[[[22], []]]
    @test merge_ds!(ds1, ds2, cotree)
    @test ds1 == [[[1], []], [[21], []], [[3], [22]]]

    # bottom cells that have the same low(e) are merged into the same cell
    cotree = Dict(11 => 1, 12 => 1, 13 => 1, 14 => 1)
    ds1 = Vector{Vector{Int}}[[[11], []]]
    ds2 = Vector{Vector{Int}}[[[12, 13], []], [[14], []]]
    @test merge_ds!(ds1, ds2, cotree)
    @test ds1 == [[[11, 12, 13], []], [[14], []]]

    # bottom cells that have the same low(e) are merged. rest ds2 edges are fused 
    cotree = Dict(11 => 1, 12 => 1, 21 => 2, 22 => 2, 3 => 3)
    ds1 = Vector{Vector{Int}}[[[11], []]]
    ds2 = Vector{Vector{Int}}[[[12], []], [[21], []], [[22], []], [[3], []]]
    @test merge_ds!(ds1, ds2, cotree)
    @test ds1 == [[[11, 12], []], [[21, 22, 3], []]]

    # ds2 edges are fused and set to the opposite side of the ds1 cell that have greater low(e)
    cotree = Dict(1 => 1, 21 => 2, 22 => 2, 23 => 2, 3 => 3)
    ds1 = Vector{Vector{Int}}[[[1], []], [[21], []], [[3], [22]]]
    ds2 = Vector{Vector{Int}}[[[23], []]]
    @test merge_ds!(ds1, ds2, cotree)
    @test ds1 == [[[1], []], [[21], []], [[3], [22, 23]]]

    # if both side of ds1 dichromatic cell have greater low(e) than ds2 bottom low(e), ds2 can not be merged
    cotree = Dict(1 => 1, 21 => 2, 22 => 2, 3 => 3, 4 => 4)
    ds1 = Vector{Vector{Int}}[[[1], []], [[21], []], [[3], [4]]]
    ds2 = Vector{Vector{Int}}[[[22], []]]
    @test !merge_ds!(ds1, ds2, cotree)

    # ds2 edges are fused and merged into top of ds1 cell if ds2 bottom has greater low(e) than that
    cotree = Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4, 5 => 5, 6 => 6)
    ds1 = Vector{Vector{Int}}[[[1], []], [[3], [2]]]
    ds2 = Vector{Vector{Int}}[[[4, 5], []], [[6], []]]
    @test merge_ds!(ds1, ds2, cotree)
    @test ds1 == [[[1], []], [[3], [2]], [[4, 5, 6], []]]
end

@testset "isplanar" begin
    @test planaritytest(complete_bipartite_graph(2,3))
    @test !planaritytest(complete_bipartite_graph(3,3))
    @test planaritytest(complete_graph(4))
    @test !planaritytest(complete_graph(5))
    @test planaritytest(ladder_graph(20))
    @test planaritytest(circular_ladder_graph(20))
    @test !planaritytest(smallgraph(:moebiuskantor))
    @test planaritytest(path_graph(20))
    disconn = SimpleGraph(Edge.([
        (1, 2), (2, 3), (3, 4), (6, 7), (6, 8),
        (6, 9), (10, 7), (10, 8), (10, 9), (11, 7), (11, 8), (11, 9)
    ]))
    @test !planaritytest(disconn)
end

@testset "isouterplanar" begin
    @test outerplanaritytest(complete_bipartite_graph(2,2))
    @test !outerplanaritytest(complete_bipartite_graph(2,3))
    @test outerplanaritytest(complete_graph(3))
    @test !outerplanaritytest(complete_graph(4))
    @test outerplanaritytest(ladder_graph(20))
    @test !outerplanaritytest(circular_ladder_graph(20))
    @test outerplanaritytest(path_graph(20))
    disconn = SimpleGraph(Edge.([
        (1, 2), (2, 3), (3, 4), (6, 7), (6, 8),
        (6, 9), (10, 7), (10, 8), (10, 9)
    ]))
    @test !outerplanaritytest(disconn)
end

end # graph.planarity
