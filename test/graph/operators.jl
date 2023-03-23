#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "operator" begin

@testset "modularproduct" begin
    g = path_graph(3)
    prod, isconn = modular_product(g, g)
    @test nv(prod) == 9
    @test ne(prod) == 10
end

@testset "line_graph" begin
    p5L, revmap, shared = line_graph(path_graph(5))
    @test nv(p5L) == 4
    @test ne(p5L) == 3
    @test length(revmap) == 4
    @test length(shared) == 3
    c5L, revmap, shared = line_graph(cycle_graph(5))
    @test nv(c5L) == 5
    @test ne(c5L) == 5
    lad5L, revmap, shared = line_graph(ladder_graph(5))
    @test nv(lad5L) == 13
    @test ne(lad5L) == 22
end

end
