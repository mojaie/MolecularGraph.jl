#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.multigraph" begin

@testset "multigraph" begin
    graph = MultiUDGraph([1,2,3,4,5], [(1,2), (3,4), (3,4), (4,5)])
    @test degree(graph, 4) == 3

end

end # graph.ugraph
