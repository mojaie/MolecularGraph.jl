#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.isomorphism.cliquebased" begin

@testset "mcis" begin
    nullg = plaingraph()
    @test length(maxcommonsubgraph(nullg, nullg, :nodeinduced)[1]) == 0

    g = plaingraph(6, [(1,2), (2,3), (2,4), (2,5), (5,6)])
    h = plaingraph(6, [(1,2), (2,3), (3,4), (3,5), (3,6)])
    @test length(maxcommonsubgraph(g, h, :nodeinduced)[1]) == 6

    g = pathgraph(4)
    h = cyclegraph(4)
    @test length(maxcommonsubgraph(g, h, :nodeinduced)[1]) == 3

    g = completegraph(7, mutable=true)
    h = unlinkedges(g, (12,))
    @test length(maxcommonsubgraph(g, h, :nodeinduced)[1]) == 6
end

@testset "mces" begin
    nullg = plaingraph()
    @test length(maxcommonsubgraph(nullg, nullg, :edgeinduced)[1]) == 0

    g = plaingraph(6, [(1,2), (2,3), (2,4), (2,5), (5,6)])
    h = plaingraph(6, [(1,2), (2,3), (3,4), (3,5), (3,6)])
    @test length(maxcommonsubgraph(g, h, :edgeinduced)[1]) == 5

    g = pathgraph(4)
    h = cyclegraph(4)
    @test length(maxcommonsubgraph(g, h, :edgeinduced)[1]) == 3

    g = completegraph(4, mutable=true)
    h = unlinkedges(g, (6,))
    @test length(maxcommonsubgraph(g, h, :edgeinduced)[1]) == 5
end


@testset "connected" begin
    g = disjointunion(completegraph(5), completegraph(4))
    h = disjointunion(completegraph(4), completegraph(3))
    @test length(maxcommonsubgraph(g, h, :nodeinduced)[1]) == 7
    @test length(maxcommonsubgraph(g, h, :nodeinduced, connected=true)[1]) == 4
end


end # graph.isomorphism.cliquebased
