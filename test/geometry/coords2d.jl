
@testset "geometry.coords2d" begin

@testset "midpoint" begin
    mp = midpoint(segment(float.([1 1; 2 2])))
    @test isapprox(rawdata(mp), [1.5, 1.5])
end

@testset "translate" begin
    seg = translate(segment(float.([1 1; 2 2])), -pi / 2, sqrt(2) / 2)
    @test isapprox(rawdata(seg), [1.5 0.5; 2.5 1.5])
end

@testset "trim" begin
    segu = trim_u(segment(float.([1 1; 5 5])), 1 / 4)
    @test isapprox(rawdata(segu), [2 2; 5 5])
    segv = trim_v(segment(float.([1 1; 5 5])), 1 / 4)
    @test isapprox(rawdata(segv), [1 1; 4 4])
    seguv = trim_uv(segment(float.([1 1; 5 5])), 1 / 4)
    @test isapprox(rawdata(seguv), [1.5 1.5; 4.5 4.5])
end

@testset "cross" begin
    @test isapprox(cross(segment([3.0 0; 0 4.0])), 12)
end

@testset "interiorangle" begin
    v = point(3.0, 0.0)
    @test isapprox(interiorangle(v, point(5.0, 0.0)), 0)
    @test isapprox(interiorangle(v, point(-20.0, 0.0)), pi)
    @test isapprox(interiorangle(v, point(3.0, -3sqrt(3))), pi / 3)
    @test isapprox(interiorangle(v, point(-3.0, 3sqrt(3))), 2pi / 3)
    @test isnan(interiorangle(v, point(0.0, 0.0)))
end

@testset "isclockwise" begin
    @test isclockwise(cyclicpath(float.([0 1; 1 0; 0 -1; -1 0])))
    @test !isclockwise(cyclicpath(float.([0 1; -1 0; 0 -1; 1 0])))
    @test isclockwise(cyclicpath(float.([1 1; 1 -1; -1 1; -1 -1]))) === nothing
    @test isclockwise(
        cyclicpath(float.([-1 1; 1 1; 1 2; 2 2; 2 -2; -2 -2; -2 2; -1 2])))
end

end # geometry.coords2d
