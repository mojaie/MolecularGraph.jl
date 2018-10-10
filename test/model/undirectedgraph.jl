
@testset "graphmodel" begin

@testset "undirectedgraph" begin
    graph = UndirectedGraph([1,2,3,4,5], [(1,2), (3,4), (4,5)])
    @test objectid(graph.nodes[3]) == objectid(graph.nodemap[3])
    @test objectid(graph.adjacency[1]) == objectid(graph.adjmap[1])
    @test objectid(graph.adjacency[5]) == objectid(graph.adjmap[5])
    node = getnode(graph, 1)
    @test typeof(node) <: Node
    @test node.index == 1
    edge = getedge(graph, 3, 4)
    @test typeof(edge) <: Edge
    @test edge.u == 3
    @test edge.v == 4
    nbr = neighbors(graph, 4)
    @test typeof(nbr) == Dict{UInt32, Edge}
    @test length(nbr) == 2
end

end # graphmodel
