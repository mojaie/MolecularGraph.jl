#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.modularproduct" begin
    g = pathgraph(3)
    h = pathgraph(3)
    prod = modularproduct(g, h)
    @test nodecount(prod) == 9
    @test edgecount(prod) == 10
end # graph.modularproduct
