
using MolecularGraph.Geometry: _point


@testset "geometry.coords3d" begin

@testset "cartesian3d" begin
    coords = cartesian3d(float.([1 2 3; 4 5 6; 7 8 9; 10 11 12]))
    @test_throws DimensionMismatch cartesian3d(float.([1 2 3 4; 5 6 7 8]))
    @test size(rawdata(coords)) == (4, 3)
    @test size(_point(coords, 3)) == (3,)
    @test sum(x_components(coords)) == 22
    @test sum(y_components(coords)) == 26
    @test sum(z_components(coords)) == 30
end

end # geometry.coords3d
