#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "traversals" begin

@testset "bfs_path" begin
    @test bfs_path([1, 1, 2, 3], 4) == [1, 2, 3, 4]
    @test bfs_path([1, 1, 2, 3], 1) == [1]
    @test bfs_path(Dict(5 => 5, 3 => 5, 7 => 3), 7) == [5, 3, 7]
end

@testset "unweighted_shortestpath" begin
    @test unweighted_shortestpath(path_graph(5), 1, 5) == [1, 2, 3, 4, 5]
    @test unweighted_shortestpath(path_graph(5), 5, 1) == [5, 4, 3, 2, 1]
    @test unweighted_shortestpath(path_graph(3), 2, 2) == [2]
    @test length(unweighted_shortestpath(cycle_graph(6), 1, 4)) == 4
    @test unweighted_shortestpath(cycle_graph(6), 1, 5) == [1, 6, 5]
    @test unweighted_shortestpath(complete_graph(5), 2, 4) == [2, 4]
    # unreachable or invalid target
    @test isempty(unweighted_shortestpath(SimpleGraph(3), 1, 3))
    @test isempty(unweighted_shortestpath(path_graph(3), 1, 5))
    g = disjoint_union(path_graph(3), path_graph(3))
    @test isempty(unweighted_shortestpath(g, 1, 6))
    @test unweighted_shortestpath(g, 4, 6) == [4, 5, 6]
end

@testset "noweight_all_distances" begin
    @test noweight_all_distances(path_graph(4), 1, 4) == [3]
    @test noweight_all_distances(cycle_graph(6), 1, 4) == [3]
    @test noweight_all_distances(cycle_graph(6), 1, 3) == [2, 4]
    @test noweight_all_distances(complete_graph(4), 1, 2) == [1, 2, 3]
    g = SimpleGraph(Edge.([(1, 2), (2, 3), (1, 4), (4, 5), (5, 3)]))
    @test noweight_all_distances(g, 1, 3) == [2, 3]
    @test isempty(noweight_all_distances(SimpleGraph(3), 1, 3))
end

end # traversals
