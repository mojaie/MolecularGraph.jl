

@testset "geometry.coords2d" begin

@testset "cartesian2d" begin
    c2d = cartesian2d(float.([1 2; 3 4; 5 6; 7 8]))
    @test size(c2d.coords) == (4, 2)
    @test length(point(c2d, 2).pos) == 2
    @test sum(x_components(c2d)) == 16
    @test sum(y_components(c2d)) == 20
end

@testset "point2d" begin
    c2d = cartesian2d(float.([1 2; 3 4; 5 6; 7 8]))
    p1 = point(c2d, 3)
    p2 = point(1.2, -4.6)
    @test x(p1) == 5
    @test y(p2) == -4.6
    @test isapprox((p1 + p2).pos, [6.2, 1.4])
    @test isapprox((p1 - p2).pos, [3.8, 10.6])
    @test isapprox((p2 * 2.0).pos, [2.4, -9.2])
    @test isapprox(cross(p1, point(c2d, 4)), -2)

    v = point(3.0, 0.0)
    @test isapprox(interiorangle(v, point(5.0, 0.0)), 0)
    @test isapprox(interiorangle(v, point(-20.0, 0.0)), pi)
    @test isapprox(interiorangle(v, point(3.0, -3sqrt(3))), pi / 3)
    @test isapprox(interiorangle(v, point(-3.0, 3sqrt(3))), 2pi / 3)
    @test isnan(interiorangle(v, point(0.0, 0.0)))
end

@testset "segment2d" begin
    c2d = cartesian2d(float.([1 2; 3 4; 5 6; 7 8]))
    seg1 = segment(c2d, 2, 4)
    seg2 = Segment2D([1.2 4.5; 3.4 7.9])
    @test size(seg2.coords) == (2, 2)
    @test u(seg1).pos == [3, 4]
    @test v(seg2).pos == [3.4, 7.9]
    @test ux(seg2) == 1.2
    @test uy(seg1) == 4
    @test vx(seg1) == 7
    @test vy(seg2) == 7.9

    @test isapprox(vector(seg2).pos, [2.2, 3.4])
    @test isapprox(midpoint(seg2).pos, [2.3, 6.2])
    tl = translate(seg1, -pi / 2, sqrt(2) / 2)
    @test isapprox(tl.coords, [3.5 3.5; 7.5 7.5])

    seg3 = Segment2D(float.([1 1; 5 5]))
    @test isapprox(trim_u(seg3, 1 / 4).coords, [2 2; 5 5])
    @test isapprox(trim_v(seg3, 1 / 4).coords, [1 1; 4 4])
    @test isapprox(trim_uv(seg3, 1 / 4).coords, [1.5 1.5; 4.5 4.5])

    @test isapprox(cross(u(seg1), v(seg1)), -4)
end

@testset "reference" begin
    c2d = cartesian2d(float.([1 2; 3 4; 5 6; 7 8]))
    seg = segment(c2d, 3, 4)
    @test seg.coords isa SubArray
    @test u(seg).pos == [5, 6]
    setcoord!(seg, 2, point(7.8, 8.3))
    @test v(seg).pos == [7.8, 8.3]
    @test point(c2d, 4).pos == [7.8, 8.3]
end

@testset "isclockwise" begin
    @test isclockwise(cartesian2d(float.([0 1; 1 0; 0 -1; -1 0])))
    @test !isclockwise(cartesian2d(float.([0 1; -1 0; 0 -1; 1 0])))
    @test isclockwise(cartesian2d(float.([1 1; 1 -1; -1 1; -1 -1]))) === nothing
    @test isclockwise(
        cartesian2d(float.([-1 1; 1 1; 1 2; 2 2; 2 -2; -2 -2; -2 2; -1 2])))
end

end # geometry.coords2d
