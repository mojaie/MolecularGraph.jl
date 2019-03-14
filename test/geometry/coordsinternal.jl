
using MolecularGraph.Geometry: _angle

@testset "geometry.coordsinternal" begin

@testset "internalcoords" begin
    emptycoords = internalcoords(10)
    @test length(emptycoords) == 10
    coords = internalcoords([
        nothing nothing nothing;
        1 nothing nothing;
        2 1 nothing;
        3 2 1
    ], [
        nothing nothing nothing;
        1.0 nothing nothing;
        1.0 2/3 nothing;
        1.0 2/3 0.0
    ])
    @test length(coords) == 4
    @test label1(coords, 2) == 1
    @test label2(coords, 3) == 1
    @test label3(coords, 4) == 1
    @test distance(coords, 2) == 1.0
    @test _angle(coords, 3) == 2 / 3
    @test dihedral(coords, 4) == 0.0
end

end # geometry.coordsinternal
