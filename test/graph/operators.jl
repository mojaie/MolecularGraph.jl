#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "operator" begin

@testset "modularproduct" begin
    nul = SimpleGraph{Int}()
    prod, isconn = modular_product(nul, nul)
    @test nv(prod) == 0
    @test ne(prod) == 0
    @test isempty(isconn)
    @test eltype(prod) === Int
    g = path_graph(3)
    prod, isconn = modular_product(g, g)
    @test nv(prod) == 9
    @test ne(prod) == 10
end

@testset "line_graph" begin
    nul, revmap, shared = line_graph(SimpleGraph{Int}())
    @test nv(nul) == 0
    @test ne(nul) == 0
    @test length(revmap) == 0
    @test length(shared) == 0
    nul, revmap, shared = line_graph(SimpleGraph(10))
    @test nv(nul) == 0
    @test ne(nul) == 0
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
