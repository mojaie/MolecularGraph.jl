#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.clique" begin

@testset "maxclique" begin
    nullg = MapUDGraph(Int[], Tuple{Int,Int}[])
    @test issetequal(maxclique(nullg), [])
    println()
    noedges = MapUDGraph(1:5, Tuple{Int,Int}[])
    @test length(maxclique(noedges)) == 1
    println()
    g1 = MapUDGraph(1:5, [(1, 2), (2, 3), (3, 4), (4, 5), (5, 1)])
    @test length(maxclique(g1)) == 2
    println()
    g2 = MapUDGraph(1:5, [(1, 2), (2, 3), (3, 1), (3, 4), (4, 5)])
    @test issetequal(maxclique(g2), [1, 2, 3])
    println()
    g3 = MapUDGraph(
        1:5, [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4), (4, 5)])
    @test issetequal(maxclique(g3), [1, 2, 3, 4])
end


end # graph.clique
