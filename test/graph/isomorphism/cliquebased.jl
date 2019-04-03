#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.isomorphism.cliquebased" begin

@testset "nodemcsclique" begin
    nullg = plaingraph()
    @test isempty(nodemcsclique(nullg, nullg))

    g = plaingraph(6, [(1,2), (2,3), (2,4), (2,5), (5,6)])
    h = plaingraph(6, [(1,2), (2,3), (3,4), (3,5), (3,6)])
    @test length(nodemcsclique(g, h)) == 6

    g = pathgraph(4)
    h = cyclegraph(4)
    @test length(nodemcsclique(g, h)) == 3

    g = completegraph(7, mutable=true)
    h = unlinkedges(g, (12,))
    @test length(nodemcsclique(g, h)) == 6
end

@testset "edgemcsclique" begin
    nullg = plaingraph()
    @test isempty(edgemcsclique(nullg, nullg))

    g = plaingraph(6, [(1,2), (2,3), (2,4), (2,5), (5,6)])
    h = plaingraph(6, [(1,2), (2,3), (3,4), (3,5), (3,6)])
    @test length(edgemcsclique(g, h)) == 5

    g = pathgraph(4)
    h = cyclegraph(4)
    @test length(edgemcsclique(g, h)) == 3

    g = completegraph(4, mutable=true)
    h = unlinkedges(g, (6,))
    @test length(edgemcsclique(g, h)) == 5
end

end # graph.isomorphism.cliquebased
