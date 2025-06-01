
@testset "util.geometry" begin

@testset "geometrybasics" begin
    # For GeometryBasics.jl interfaces compatiblity
    # Points
    p1 = Point(5, 6)
    p2 = Point(1.2, -4.6)
    @test p1[1] == 5
    @test p2[2] == -4.6
    @test isapprox(p1 + p2, Point(6.2, 1.4))
    @test isapprox(p1 - p2, Point(3.8, 10.6))
    @test isapprox(p2 * 2.0, p2 / 0.5)
    @test isapprox(cross(p1, Point(7, 8)), -2)

    # Lines
    seg1 = Line(Point(3, 4), Point(7, 8))
    seg2 = Line(Point(1.2, 4.5), Point(3.4, 7.9))
    @test seg1[1] == Point(3.0, 4.0)
    @test seg2[2] == Point(3.4, 7.9)
    @test isapprox(seg2[1] + (seg2[2] - seg2[1]) / 2, Point(2.3, 6.2))
end

@testset "geometry" begin
    v = Point(3.0, 0.0)
    @test isapprox(interiorangle(v, Point(5.0, 0.0)), 0)
    @test isapprox(interiorangle(v, Point(-20.0, 0.0)), pi)
    @test isapprox(interiorangle(v, Point(3.0, -3sqrt(3))), pi / 3)
    @test isapprox(interiorangle(v, Point(-3.0, 3sqrt(3))), 2pi / 3)
    @test isnan(interiorangle(v, Point(0.0, 0.0)))

    seg1 = Line(Point(3, 4), Point(7, 8))
    tl = translate(seg1, -pi / 2, sqrt(2) / 2)
    @test isapprox(tl[1], Point(3.5, 3.5))
    @test isapprox(tl[2], Point(7.5, 7.5))

    seg3 = Line(Point(1, 1), Point(5, 5))
    tu = trim_u(seg3, sqrt(2))
    @test isapprox(tu[1], Point(2, 2))
    @test isapprox(tu[2], Point(5, 5))
    tv = trim_v(seg3, sqrt(2))
    @test isapprox(tv[1], Point(1, 1))
    @test isapprox(tv[2], Point(4, 4))
    tuv = trim_uv(seg3, sqrt(2) / 2)
    @test isapprox(tuv[1], Point(1.5, 1.5))
    @test isapprox(tuv[2], Point(4.5, 4.5))

    @test isclockwise(Point2d.([[0, 1], [1, 0], [0, -1], [-1, 0]]))
    @test !isclockwise(Point2d.([[0, 1], [-1, 0], [0, -1], [1, 0]]))
    @test isclockwise(Point2d.([[1, 1], [1, -1], [-1, 1], [-1, -1]])) === nothing
    @test isclockwise(Point2d.([[-1, 1], [1, 1], [1, 2], [2, 2], [2, -2], [-2, -2], [-2, 2], [-1, 2]]))
end


end # util.geometry
