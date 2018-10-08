@testset "undirectedgraph" begin
    graph = UndirectedGraph([1,2,3,4,5], [(1,2), (3,4), (4,5)])
    node = getnode(graph, 1)
    @test typeof(node) <: Node
    @test node.index == 1
    edge = getedge(graph, 3, 4)
    @test typeof(edge) <: Edge
    @test edge.u == 3
    @test edge.v == 4
    nbr = getneighbors(graph, 4)
    @test typeof(nbr) == Dict{UInt32, Edge}
    @test length(nbr) == 2
end
