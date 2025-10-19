#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


@testset "coords.descriptor" begin
    desc = MolDescriptor{Int}()
    push!(desc.draw2d_bond_style, Draw2dBondStyle([:revup, :none, :cis_trans, :none]))
    remap!(
        Val(:draw2d_bond_style), desc, [1, 2, 5, 4],
        Edge.([(1, 2), (1, 3), (1, 4), (3, 5)])
    )
    @test desc.draw2d_bond_style[1] == [:revup, :cis_trans]

    desc = MolDescriptor{Int}()
    push!(desc.draw2d_bond_style, Draw2dBondStyle([:none, :none, :none, :unspecified, :up, :down]))
    remap!(
        Val(:draw2d_bond_style), desc, [1, 7, 3, 6, 5],
        Edge.([(1, 2), (2, 3), (3, 4), (3, 5), (3, 6), (5, 7)])
    )
    @test desc.draw2d_bond_style[1] == Draw2dBondStyle([:revdown, :up, :unspecified])
    dump = JSON.json(desc.draw2d_bond_style)
    @test JSON.parse(dump, Vector{Draw2dBondStyle}) == desc.draw2d_bond_style
end


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

    # TODO: coordgen error : CHMMol-> sketcher stereo error: wrong number for RSpriorities
    # Takes about 4 sec to calculate
    # thiostrepton = smilestomol(raw"C[C@@H](O)[C@@](C)(O)[C@@H](C2=NC6=CS2)NC(C1CSC(/C(NC([C@@H](NC(C3=CS[C@]([C@]4(NC([C@@H](NC(C(NC([C@H](C)NC%10=O)=O)=C)=O)C)=O)C(C5=CSC([C@H]([C@@H](C)OC(C8=CC([C@H](C)O)=C(C=CC(N[C@]([H])%10[C@H](CC)C)[C@@H]9O)C9=N8)=O)NC6=O)=N5)N=C(C7=NC(C(NC(C(NC(C(N)=O)=C)=O)=C)=O)=CS7)CC4)=N3)=O)[C@@H](C)O)=O)=C/C)=N1)=O")
    # @test get_prop(thiostrepton, :stereobond)[Edge(21 => 123)] == (20, 124, false)
end
