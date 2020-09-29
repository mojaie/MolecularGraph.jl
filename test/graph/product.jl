#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.product" begin

    @testset "modularproduct" begin
        g = pathgraph(3)
        prod = modularproduct(g, g)
        @test nodecount(prod) == 9
        @test edgecount(prod) == 10
    end

    @testset "cartesian.product" begin
        g = pathgraph(5)
        prod = cartesianproduct(g, g)
        grid = squaregrid(5, 5)
        @test isexactmatch(prod, grid)
    end
end # graph.product
