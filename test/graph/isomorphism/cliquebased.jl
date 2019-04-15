#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.isomorphism.cliquebased" begin

@testset "mcis" begin
    nullg = plaingraph()
    @test mcissize(nullg, nullg) == 0

    g = plaingraph(6, [(1,2), (2,3), (2,4), (2,5), (5,6)])
    h = plaingraph(6, [(1,2), (2,3), (3,4), (3,5), (3,6)])
    @test mcissize(g, h) == 6

    g = pathgraph(4)
    h = cyclegraph(4)
    @test mcissize(g, h) == 3

    g = completegraph(7, mutable=true)
    h = unlinkedges(g, (12,))
    @test mcissize(g, h) == 6
end

@testset "mces" begin
    nullg = plaingraph()
    @test mcessize(nullg, nullg) == 0

    g = plaingraph(6, [(1,2), (2,3), (2,4), (2,5), (5,6)])
    h = plaingraph(6, [(1,2), (2,3), (3,4), (3,5), (3,6)])
    @test mcessize(g, h) == 5

    g = pathgraph(4)
    h = cyclegraph(4)
    @test mcessize(g, h) == 3

    g = completegraph(4, mutable=true)
    h = unlinkedges(g, (6,))
    @test mcessize(g, h) == 5
end

end # graph.isomorphism.cliquebased
