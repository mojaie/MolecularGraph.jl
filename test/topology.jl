
@testset "topology.Ring" begin
    a = GraphMol.Ring([1, 2, 3, 4, 5])
    @test typeof(a) == GraphMol.Ring
    b = GraphMol.Ring([5, 1, 2, 3, 4])
    @test a == b
    b = GraphMol.Ring([1, 2, 3, 4, 5, 6])
    @test a != b
    b = GraphMol.Ring([1, 2, 3, 5, 4])
    @test a != b
    r = GraphMol.Ring([5, 3, 4, 2, 1])
    @test r.arr == UInt16[1, 2, 4, 3, 5]
    a = GraphMol.Ring([5, 6, 7, 8, 9, 1])
    @test a.arr == UInt16[1, 5, 6, 7, 8, 9]
end


@testset "topology.resolve_inclusion" begin
    a = GraphMol.Ring([1, 2, 3, 4, 5])
    b = GraphMol.Ring([1, 2, 3, 4, 5, 6, 7, 8, 9])
    @test GraphMol.resolve_inclusion(a, b)[2].arr == UInt16[1, 5, 6, 7, 8, 9]
    b = GraphMol.Ring([1, 6, 7, 8, 9])
    @test GraphMol.resolve_inclusion(a, b) == nothing
    b = GraphMol.Ring([3, 4, 5, 6, 7, 8, 9])
    @test GraphMol.resolve_inclusion(a, b) == nothing
    b = GraphMol.Ring([2, 3, 4, 5, 6, 7, 8])
    @test GraphMol.resolve_inclusion(a, b)[2].arr == UInt16[1, 2, 8, 7, 6, 5]
    a = GraphMol.Ring([3, 4, 5, 6, 7, 8, 9, 10])
    b = GraphMol.Ring([1, 2, 8, 7, 6, 5, 4])
    @test GraphMol.resolve_inclusion(a, b)[1].arr == UInt16[1, 2, 8, 9, 10, 3, 4]
    a = GraphMol.Ring([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
    b = GraphMol.Ring([8, 9, 10, 11, 12])
    @test GraphMol.resolve_inclusion(a, b)[1].arr == [1, 2, 3, 4, 5, 6, 7, 8, 12, 11]
end


@testset "topology.topology" begin
    RESOURCE_DIR = joinpath(dirname(@__FILE__), "..", "_resources", "DrugBank")
    phe = GraphMol.loadsdfmol(open(joinpath(RESOURCE_DIR, "Phe.mol")))
end
