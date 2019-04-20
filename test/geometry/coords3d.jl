

@testset "geometry.coords3d" begin

@testset "cartesian3d" begin
    c3d = cartesian3d(float.([1 2 3; 4 5 6; 7 8 9; 10 11 12]))
    @test size(c3d.coords) == (4, 3)
    @test sum(x_components(c3d)) == 22
    @test sum(y_components(c3d)) == 26
    @test sum(z_components(c3d)) == 30
end

end # geometry.coords3d
