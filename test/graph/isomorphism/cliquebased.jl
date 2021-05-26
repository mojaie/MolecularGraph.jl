#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.isomorphism.cliquebased" begin

@testset "mcis" begin
    nullg = plaingraph()
    @test size(maxcommonsubgraph(nullg, nullg, :nodeinduced)) == 0

    g = plaingraph(6, [(1,2), (2,3), (2,4), (2,5), (5,6)])
    h = plaingraph(6, [(1,2), (2,3), (3,4), (3,5), (3,6)])
    @test size(maxcommonsubgraph(g, h, :nodeinduced)) == 6

    g = pathgraph(4)
    h = cyclegraph(4)
    @test size(maxcommonsubgraph(g, h, :nodeinduced)) == 3

    g = completegraph(7, mutable=true)
    h = unlinkedges(g, (12,))
    @test size(maxcommonsubgraph(g, h, :nodeinduced)) == 6
end

@testset "mces" begin
    nullg = plaingraph()
    @test size(maxcommonsubgraph(nullg, nullg, :edgeinduced)) == 0

    g = plaingraph(6, [(1,2), (2,3), (2,4), (2,5), (5,6)])
    h = plaingraph(6, [(1,2), (2,3), (3,4), (3,5), (3,6)])
    @test size(maxcommonsubgraph(g, h, :edgeinduced)) == 5

    g = pathgraph(4)
    h = cyclegraph(4)
    @test size(maxcommonsubgraph(g, h, :edgeinduced)) == 3

    g = completegraph(4, mutable=true)
    h = unlinkedges(g, (6,))
    @test size(maxcommonsubgraph(g, h, :edgeinduced)) == 5
end


@testset "connected" begin
    g = disjointunion(completegraph(5), completegraph(4))
    h = disjointunion(completegraph(4), completegraph(3))
    @test size(maxcommonsubgraph(g, h, :nodeinduced)) == 7
    @test size(maxcommonsubgraph(g, h, :nodeinduced, connected=true)) == 4
end


end # graph.isomorphism.cliquebased
