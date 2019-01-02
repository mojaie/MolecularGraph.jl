#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.linegraph" begin

@testset "linegraph" begin
    G = VectorUDGraph(6, [(1,2), (2,3), (2,4), (2,5), (5,6)])
    L = linegraph(G)
    @test nodecount(L) == 5
    @test edgecount(L) == 7
end

end # graph.linegraph
