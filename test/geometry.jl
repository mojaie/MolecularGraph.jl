@testset "geometry" begin
    @test distance((2.0, 2.0), (5.0, 6.0)) ≈ 5
    @test distance([2, 2], [2, 2]) ≈ 0
    @test rotate((4.0, 4.0), pi / 2) ≈ Point2D(-4, 4)
    @test cross2d((3.0, 0), (0, 4.0)) ≈ 12
    @test interiorangle((3.0, 0), (5.0, 0)) ≈ 0
    @test interiorangle((3.0, 0), (-20.0, 0)) ≈ pi
    @test interiorangle((3.0, 0), (3.0, -3sqrt(3))) ≈ pi / 3
    @test interiorangle((3.0, 0), (-3.0, 3sqrt(3))) ≈ 2pi / 3
    @test isnan(interiorangle((3.0, 0), (0, 0)))
    @test isclockwise([(0, 1), (1, 0), (0, -1), (-1, 0)])
    @test !isclockwise([(0, 1), (-1, 0), (0, -1), (1, 0)])
    @test isnan(isclockwise([(1, 1), (1, -1), (-1, 1), (-1, -1)]))
    @test isclockwise([
        (-1, 1), (1, 1), (1, 2), (2, 2), (2, -2), (-2, -2), (-2, 2), (-1, 2)
    ])
    pm1 = parallelmove((1, 1), (2, 2), -pi / 2, sqrt(2) / 2)
    @test pm1[1] ≈ Point2D(1.5, 0.5)
    @test pm1[2] ≈ Point2D(2.5, 1.5)
    pt1 = paralleltrim((1, 1), (5, 5), 1 / 4, 1)
    @test pt1[1] ≈ Point2D(1, 1)
    @test pt1[2] ≈ Point2D(4, 4)
    pt2 = paralleltrim((1, 1), (5, 5), 1 / 4, 2)
    @test pt2[1] ≈ Point2D(2, 2)
    @test pt2[2] ≈ Point2D(5, 5)
    pt3 = paralleltrim((1, 1), (5, 5), 1 / 4)
    @test pt3[1] ≈ Point2D(1.5, 1.5)
    @test pt3[2] ≈ Point2D(4.5, 4.5)
end
