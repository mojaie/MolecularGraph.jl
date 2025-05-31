
@testset "geometry.embed2d" begin

@testset "embed2d" begin
    return
    mol = smilestomol("CCC1CCCC1C")
    zmatrix = coords2d(mol)
    display(zmatrix)
    println("")
    coords3d = cartesian(zmatrix)
    display(coords3d)
    println("")
    for (i, atom) in enumerate(atomvector(mol))
        atom.coords = tuple(coords3d[i, 1:2]...)
    end
    ASSETS_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets")
    dest = open(joinpath(ASSETS_DIR, "image", "test.svg"), "w")
    write(dest, drawsvg(mol, 200, 200))
end

end # geometry.embed2d
