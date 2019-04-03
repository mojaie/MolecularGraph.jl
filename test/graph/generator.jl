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
    @test distance(P100, 1, 100) == 99
end

@testset "cyclegraph" begin
    @test_throws DomainError cyclegraph(2)
    C3 = cyclegraph(3)
    @test nodecount(C3) == edgecount(C3)
    C100 = cyclegraph(100)
    @test nodecount(C100) == edgecount(C100)
    @test distance(C100, 1, 100) == 1
    @test distance(C100, 1, 50) == 49
end

@testset "completebipartite" begin
    @test_throws DomainError completebipartite(0, 1)
    @test_throws DomainError completebipartite(1, 0)
    K1_1 = completebipartite(1, 1)
    @test nodecount(K1_1) == 2
    @test edgecount(K1_1) == 1
    K3_3 = completebipartite(3, 3)
    @test nodecount(K3_3) == 6
    @test edgecount(K3_3) == 9
end

@testset "completegraph" begin
    @test_throws DomainError completegraph(-1)
    K0 = completegraph(0)
    @test nodecount(K0) == 0
    @test edgecount(K0) == 0
    K1 = completegraph(1)
    @test nodecount(K1) == 1
    @test edgecount(K1) == 0
    K2 = completegraph(2)
    @test nodecount(K2) == 2
    @test edgecount(K2) == 1
    K3 = completegraph(3)
    @test nodecount(K3) == 3
    @test edgecount(K3) == 3
    K5 = completegraph(5)
    @test nodecount(K5) == 5
    @test edgecount(K5) == 10
    K20 = completegraph(20)
    @test nodecount(K20) == 20
    @test edgecount(K20) == 190
end

@testset "laddergraph" begin
    @test_throws DomainError laddergraph(0)
    L3 = laddergraph(3)
    @test edgecount(L3) == 7 # 3n - 2
    L100 = laddergraph(100)
    @test edgecount(L100) == 298 # 3n - 2
    @test distance(L100, 1, 100) == 50
end

@testset "circularladder" begin
    @test_throws DomainError circularladder(2)
    CL3 = circularladder(3)
    @test edgecount(CL3) == 9 # 3n
    CL100 = circularladder(100)
    @test edgecount(CL100) == 300 # 3n
    @test distance(CL100, 1, 50) == 25
    @test distance(CL100, 1, 200) == 2
end

@testset "moebiusladder" begin
    @test_throws DomainError moebiusladder(2)
    ML3 = moebiusladder(3)
    @test edgecount(ML3) == 9 # 3n
    ML100 = moebiusladder(100)
    @test edgecount(ML100) == 300 # 3n
    @test distance(ML100, 1, 50) == 25
    @test distance(ML100, 1, 200) == 1
end

end # graph.generator
