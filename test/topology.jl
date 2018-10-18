
@testset "topology" begin

@testset "Ring" begin
    a = Ring([1, 2, 3, 4, 5])
    @test typeof(a) == Ring
    b = Ring([5, 1, 2, 3, 4])
    @test a == b
    b = Ring([1, 2, 3, 4, 5, 6])
    @test a != b
    b = Ring([1, 2, 3, 5, 4])
    @test a != b
    r = Ring([5, 3, 4, 2, 1])
    @test r.arr == UInt16[1, 2, 4, 3, 5]
    a = Ring([5, 6, 7, 8, 9, 1])
    @test a.arr == UInt16[1, 5, 6, 7, 8, 9]
end


@testset "resolve_inclusion" begin
    a = Ring([1, 2, 3, 4, 5])
    b = Ring([1, 2, 3, 4, 5, 6, 7, 8, 9])
    @test resolve_inclusion(a, b)[2].arr == UInt16[1, 5, 6, 7, 8, 9]
    b = Ring([1, 6, 7, 8, 9])
    @test resolve_inclusion(a, b) == nothing
    b = Ring([3, 4, 5, 6, 7, 8, 9])
    @test resolve_inclusion(a, b) == nothing
    b = Ring([2, 3, 4, 5, 6, 7, 8])
    @test resolve_inclusion(a, b)[2].arr == UInt16[1, 2, 8, 7, 6, 5]
    a = Ring([3, 4, 5, 6, 7, 8, 9, 10])
    b = Ring([1, 2, 8, 7, 6, 5, 4])
    @test resolve_inclusion(a, b)[1].arr == UInt16[1, 2, 8, 9, 10, 3, 4]
    a = Ring([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
    b = Ring([8, 9, 10, 11, 12])
    @test resolve_inclusion(a, b)[1].arr == [1, 2, 3, 4, 5, 6, 7, 8, 12, 11]
end

"""
@testset "topology" begin
    RESOURCE_DIR = joinpath(dirname(@__FILE__), "..", "_resources", "DrugBank")
    phe = loadsdfmol(open(joinpath(RESOURCE_DIR, "Phe.mol")))
    @test phe.rings[1].arr == UInt16[6, 7, 10, 12, 11, 8]
    @test phe.scaffolds == [[1]]
    @test length(phe.isolated) == 0
    prem = loadsdfmol(open(joinpath(RESOURCE_DIR, "Premarin.mol")))
    @test prem.scaffolds == [[1, 2, 3, 4]]
    @test prem.isolated == [[28]]
end
"""

end # topology
