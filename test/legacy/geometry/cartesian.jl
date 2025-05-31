

@testset "geometry.cartesian" begin

@testset "point2d" begin
    coords = float.([1 2; 3 4; 5 6; 7 8])
    p1 = Point2D(coords, 3)
    p2 = Point2D(1.2, -4.6)
    @test p1.x == 5
    @test p2.y == -4.6
    @test isapprox(toarray(p1 + p2), [6.2, 1.4])
    @test isapprox(toarray(p1 - p2), [3.8, 10.6])
    @test isapprox(toarray(p2 * 2.0), [2.4, -9.2])
    @test isapprox(cross2d(p1, Point2D(coords, 4)), -2)

    v = Point2D(3.0, 0.0)
    @test isapprox(interiorangle(v, Point2D(5.0, 0.0)), 0)
    @test isapprox(interiorangle(v, Point2D(-20.0, 0.0)), pi)
    @test isapprox(interiorangle(v, Point2D(3.0, -3sqrt(3))), pi / 3)
    @test isapprox(interiorangle(v, Point2D(-3.0, 3sqrt(3))), 2pi / 3)
    @test isnan(interiorangle(v, Point2D(0.0, 0.0)))
end

@testset "segment2d" begin
    coords = float.([1 2; 3 4; 5 6; 7 8])
    seg1 = Segment{Point2D}(coords, 2, 4)
    seg2 = Segment(Point2D(1.2, 4.5), Point2D(3.4, 7.9))
    @test toarray(seg1.u) == [3, 4]
    @test toarray(seg2.v) == [3.4, 7.9]
    @test seg2.u.x == 1.2
    @test seg1.u.y == 4
    @test seg1.v.x == 7
    @test seg2.v.y == 7.9

    @test isapprox(toarray(seg2.v - seg2.u), [2.2, 3.4])
    @test isapprox(toarray(midpoint(seg2)), [2.3, 6.2])
    tl = translate(seg1, -pi / 2, sqrt(2) / 2)
    @test isapprox(toarray(tl), [3.5 3.5; 7.5 7.5])

    seg3 = Segment(Point2D(1, 1), Point2D(5, 5))
    @test isapprox(toarray(trim_u(seg3, 1 / 4)), [2 2; 5 5])
    @test isapprox(toarray(trim_v(seg3, 1 / 4)), [1 1; 4 4])
    @test isapprox(toarray(trim_uv(seg3, 1 / 4)), [1.5 1.5; 4.5 4.5])

    @test isapprox(cross2d(seg1.u, seg1.v), -4)
end

@testset "cartesian2d" begin
    coords = float.([1 2; 3 4; 5 6; 7 8])
    @test sum(x_components(coords)) == 16
    @test sum(y_components(coords)) == 20
    seg = toarray(coords, 3, 4)
    @test seg isa SubArray
end

@testset "cartesian3d" begin
    coords = float.([1 2 3; 4 5 6; 7 8 9; 10 11 12])
    @test sum(x_components(coords)) == 22
    @test sum(y_components(coords)) == 26
    @test sum(z_components(coords)) == 30
end

@testset "isclockwise" begin
    @test isclockwise(float.([0 1; 1 0; 0 -1; -1 0]))
    @test !isclockwise(float.([0 1; -1 0; 0 -1; 1 0]))
    @test isclockwise(float.([1 1; 1 -1; -1 1; -1 -1])) === nothing
    @test isclockwise(float.([-1 1; 1 1; 1 2; 2 2; 2 -2; -2 -2; -2 2; -1 2]))
end

end # geometry.cartesian
