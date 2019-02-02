#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.generator" begin

@testset "pathgraph" begin
    @test_throws DomainError pathgraph(1)
    P2 = pathgraph(2)
    @test nodecount(P2) == edgecount(P2) + 1
    P100 = pathgraph(100)
    @test nodecount(P100) == edgecount(P100) + 1
    @test length(shortestpath(P100, 1, 100)) == 100
end

@testset "cyclegraph" begin
    @test_throws DomainError cyclegraph(2)
    C3 = cyclegraph(3)
    @test nodecount(C3) == edgecount(C3)
    C100 = cyclegraph(100)
    @test nodecount(C100) == edgecount(C100)
    @test length(shortestpath(C100, 1, 100)) == 2
    @test length(shortestpath(C100, 1, 50)) == 50
end

@testset "completegraph" begin
    @test_throws DomainError completegraph(-1)
    cg0 = completegraph(0)
    @test nodecount(cg0) == 0
    @test edgecount(cg0) == 0
    cg1 = completegraph(1)
    @test nodecount(cg1) == 1
    @test edgecount(cg1) == 0
    cg2 = completegraph(2)
    @test nodecount(cg2) == 2
    @test edgecount(cg2) == 1
    cg3 = completegraph(3)
    @test nodecount(cg3) == 3
    @test edgecount(cg3) == 3
    cg5 = completegraph(5)
    @test nodecount(cg5) == 5
    @test edgecount(cg5) == 10
    cg20 = completegraph(20)
    @test nodecount(cg20) == 20
    @test edgecount(cg20) == 190
end

end # graph.generator
