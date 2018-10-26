
@testset "geometry" begin

@testset "cross" begin
    @test cross(SVector(3.0, 0), SVector(0, 4.0)) ≈ 12
end

@testset "interiorangle" begin
    v = SVector(3.0, 0)
    @test interiorangle(v, SVector(5.0, 0)) ≈ 0
    @test interiorangle(v, SVector(-20.0, 0)) ≈ pi
    @test interiorangle(v, SVector(3.0, -3sqrt(3))) ≈ pi / 3
    @test interiorangle(v, SVector(-3.0, 3sqrt(3))) ≈ 2pi / 3
    @test isnan(interiorangle(v, SVector(0, 0)))
end

@testset "isclockwise" begin
    @test isclockwise(SVector{2}[[0, 1], [1, 0], [0, -1], [-1, 0]])
    @test !isclockwise(SVector{2}[[0, 1], [-1, 0], [0, -1], [1, 0]])
    @test isnan(isclockwise(SVector{2}[[1, 1], [1, -1], [-1, 1], [-1, -1]]))
    @test isclockwise(SVector{2}[
        [-1, 1], [1, 1], [1, 2], [2, 2], [2, -2], [-2, -2], [-2, 2], [-1, 2]
    ])
end

@testset "translate" begin
    uv = @SMatrix [1 1; 2 2]
    uv1 = translate(uv, -pi / 2, sqrt(2) / 2)
    @test uv1[1, :] ≈ [1.5, 0.5]
    @test uv1[2, :] ≈ [2.5, 1.5]
end

@testset "trim" begin
    uv = @SMatrix [1 1; 5 5]
    pt1 = trimV(uv, 1 / 4)
    @test pt1[1, :] ≈ [1, 1]
    @test pt1[2, :] ≈ [4, 4]
    pt2 = trimU(uv, 1 / 4)
    @test pt2[1, :] ≈ [2, 2]
    @test pt2[2, :] ≈ [5, 5]
    pt3 = trimUV(uv, 1 / 4)
    @test pt3[1, :] ≈ [1.5, 1.5]
    @test pt3[2, :] ≈ [4.5, 4.5]
end

end # geometry
