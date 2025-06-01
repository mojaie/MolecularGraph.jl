#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "coords" begin
    # Rifampicin
    mol = smilestomol(raw"CN1CCN(CC1)/N=C/c2c(O)c3c5C(=O)[C@@]4(C)O/C=C/[C@H](OC)[C@@H](C)[C@@H](OC(C)=O)[C@H](C)[C@H](O)[C@H](C)[C@@H](O)[C@@H](C)\C=C\C=C(\C)C(=O)Nc2c(O)c3c(O)c(C)c5O4")
    stereocenter_from_smiles!(mol)
    stereobond_from_smiles!(mol)
    coords, styles = coordgen(mol)
    # just check size
    @test length(coords) == nv(mol)
    @test length(styles) == ne(mol)
    # println(coords)
    # println(styles)

    # maybe bug in coordgenlib, wedge isReversed may not work
    # jsm = smilestomol(raw"CC/C=C\C[C@@H]1[C@H](CCC1=O)CC(=O)O")
end
