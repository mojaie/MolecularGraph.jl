#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.multigraph" begin

@testset "multigraph" begin
    graph = multigraph(5, [(1,2), (3,4), (3,4), (4,5)])
    @test degree(graph, 4) == 3
    @test length(adjacencies(graph, 4)) == 2
    @test length(incidences(graph, 4)) == 3
end

end # graph.ugraph
