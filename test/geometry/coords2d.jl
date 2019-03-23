
using MolecularGraph.MolecularGraphGeometry:
    _point, _u, _v, _vector, _midpoint, _translate, _trim_u, _trim_v, _trim_uv

@testset "geometry.coords2d" begin

@testset "cartesian2d" begin
    coords = cartesian2d(float.([1 2; 3 4; 5 6; 7 8]))
    @test_throws DimensionMismatch cartesian2d(float.([1 2 3; 4 5 6]))
    @test size(rawdata(coords)) == (4, 2)
    @test size(_point(coords, 2)) == (2,)
    @test sum(x_components(coords)) == 16
    @test sum(y_components(coords)) == 20
end

@testset "point2d" begin
    coords = cartesian2d(float.([1 2; 3 4; 5 6; 7 8]))
    p1 = point(coords, 3)
    p2 = point(1.2, -4.6)
    @test size(rawdata(p2)) == (2,)
    @test y(p2) == -4.6
    @test x(p1) == 5

    (x1, y1) = Formatting.fmt("05.2f", p2)
    @test x1 == "01.20"
    @test y1 == "-4.60"

    @test isapprox(cross(p1, point(coords, 4)), -2)

    v = point(3.0, 0.0)
    @test isapprox(interiorangle(v, point(5.0, 0.0)), 0)
    @test isapprox(interiorangle(v, point(-20.0, 0.0)), pi)
    @test isapprox(interiorangle(v, point(3.0, -3sqrt(3))), pi / 3)
    @test isapprox(interiorangle(v, point(-3.0, 3sqrt(3))), 2pi / 3)
    @test isnan(interiorangle(v, point(0.0, 0.0)))
end

@testset "segment2d" begin
    coords = cartesian2d(float.([1 2; 3 4; 5 6; 7 8]))
    seg1 = segment(coords, 2, 4)
    seg2 = segment([1.2 4.5; 3.4 7.9])
    @test size(rawdata(seg2)) == (2, 2)
    @test _u(seg1) == [3, 4]
    @test _v(seg2) == [3.4, 7.9]
    @test ux(seg2) == 1.2
    @test uy(seg1) == 4
    @test vx(seg1) == 7
    @test vy(seg2) == 7.9

    ((x1, y1), (x2, y2)) = Formatting.fmt("05.2f", seg2)
    @test y1 == "04.50"
    @test x2 == "03.40"
    @test isapprox(_vector(seg2), [2.2, 3.4])
    @test isapprox(_midpoint(seg2), [2.3, 6.2])
    @test isapprox(_translate(seg1, -pi / 2, sqrt(2) / 2), [3.5 3.5; 7.5 7.5])

    seg3 = segment(float.([1 1; 5 5]))
    @test isapprox(_trim_u(seg3, 1 / 4), [2 2; 5 5])
    @test isapprox(_trim_v(seg3, 1 / 4), [1 1; 4 4])
    @test isapprox(_trim_uv(seg3, 1 / 4), [1.5 1.5; 4.5 4.5])

    @test isapprox(cross(seg1), -4)
end

@testset "reference" begin
    coords = cartesian2d(float.([1 2; 3 4; 5 6; 7 8]))
    seg = segment(coords, 3, 4)
    pos = point(seg, 1)
    @test rawdata(pos) == [5, 6]
    setcoord!(seg, point(7.8, 8.3), 2)
    @test _point(seg, 2) == [7.8, 8.3]
    @test _point(coords, 4) == [7, 8]
end

@testset "isclockwise" begin
    @test isclockwise(cyclicpath(float.([0 1; 1 0; 0 -1; -1 0])))
    @test !isclockwise(cyclicpath(float.([0 1; -1 0; 0 -1; 1 0])))
    @test isclockwise(cyclicpath(float.([1 1; 1 -1; -1 1; -1 -1]))) === nothing
    @test isclockwise(
        cyclicpath(float.([-1 1; 1 1; 1 2; 2 2; 2 -2; -2 -2; -2 2; -1 2])))
end

end # geometry.coords2d
