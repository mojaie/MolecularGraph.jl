
@testset "geometry" begin

@testset "distance" begin
    @test distance((2.0, 2.0), (5.0, 6.0)) ≈ 5
    @test distance([2, 2], [2, 2]) ≈ 0
end

@testset "rotate" begin
    @test rotate((4.0, 4.0), pi / 2) ≈ Point2D(-4, 4)
end

@testset "cross2d" begin
    @test cross2d((3.0, 0), (0, 4.0)) ≈ 12
end

@testset "interiorangle" begin
    @test interiorangle((3.0, 0), (5.0, 0)) ≈ 0
    @test interiorangle((3.0, 0), (-20.0, 0)) ≈ pi
    @test interiorangle((3.0, 0), (3.0, -3sqrt(3))) ≈ pi / 3
    @test interiorangle((3.0, 0), (-3.0, 3sqrt(3))) ≈ 2pi / 3
    @test isnan(interiorangle((3.0, 0), (0, 0)))
end

@testset "isclockwise" begin
    @test isclockwise([(0, 1), (1, 0), (0, -1), (-1, 0)])
    @test !isclockwise([(0, 1), (-1, 0), (0, -1), (1, 0)])
    @test isnan(isclockwise([(1, 1), (1, -1), (-1, 1), (-1, -1)]))
    @test isclockwise([
        (-1, 1), (1, 1), (1, 2), (2, 2), (2, -2), (-2, -2), (-2, 2), (-1, 2)
    ])
end

@testset "translate" begin
    seg = segment((1, 1), (2, 2))
    pm1 = translate(seg, -pi / 2, sqrt(2) / 2)
    @test pm1.u ≈ Point2D(1.5, 0.5)
    @test pm1.v ≈ Point2D(2.5, 1.5)
end

@testset "trim" begin
    pt1 = trim_v(segment((1, 1), (5, 5)), 1 / 4)
    @test pt1.u ≈ Point2D(1, 1)
    @test pt1.v ≈ Point2D(4, 4)
    pt2 = trim_u(segment((1, 1), (5, 5)), 1 / 4)
    @test pt2.u ≈ Point2D(2, 2)
    @test pt2.v ≈ Point2D(5, 5)
    pt3 = trim_uv(segment((1, 1), (5, 5)), 1 / 4)
    @test pt3.u ≈ Point2D(1.5, 1.5)
    @test pt3.v ≈ Point2D(4.5, 4.5)
end

end # geometry
