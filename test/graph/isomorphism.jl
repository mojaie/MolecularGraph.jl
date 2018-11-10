
@testset "isomorphism" begin

@testset "subgraph" begin
    g = UDGraph(6, [(1,2), (2,3), (2,4), (2,5), (5,6)])
    h = UDGraph(6, [(1,2), (2,3), (3,4), (3,5), (3,6)])
    @test subgraph_is_isomorphic(g, h)
    @test is_isomorphic(g, h)
    h2 = UDGraph(7, [(1,2), (2,3), (3,4), (3,5), (3,6), (6,7)])
    @test subgraph_is_isomorphic(h2, g)
    @test !is_isomorphic(h2, g)
    h3 = UDGraph(6, [(1,2), (2,3), (3,4), (3,5), (3,6), (6,1)])
    @test !subgraph_is_isomorphic(h3, g)
end

end # graphmodel
