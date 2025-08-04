
@testset "virtualatom" begin

@testset "virtualatom" begin
    mol = smilestomol(GeneralMolGraph, "c1ccccc1C")
    set_prop!(mol, 7, VirtualAtom("R"))
    @test valence(mol)[7] == 1
    @test atom_symbol(mol.vprops[7]) === :R
    @test molecular_formula(mol) == "C6H5R"
    @test is_ring_aromatic(mol)[1]
    set_prop!(mol, 7, VirtualAtom("<sub>R</sub>"))
    @test atom_symbol(mol.vprops[7]) === Symbol("&lt;sub&gt;R&lt;/sub&gt;")

    mol = smilestomol(GeneralMolGraph, "n1cncc1"; on_update=default_on_update!, initialized=true)
    set_prop!(mol, 1, HydrogenatedAtom(SMILESAtom(:N), 1))
    set_descriptor!(mol, :sssr, mincyclebasis(mol.graph))  # initialize SSSR
    smiles_on_update!(mol)
    @test total_hydrogens(mol) == [1, 1, 0, 1, 1]
    @test is_ring_aromatic(mol)[1]

    mol = sdftomol(GeneralMolGraph, joinpath(dirname(@__FILE__), "../assets/test/demo.mol"))
    set_prop!(mol, 6, FormulaGroup(
        Dict(:H => 7, :C => 2, :N => 1), 1, Vector{Tuple{Symbol, String}}[]))
    @test valence(mol)[6] == 1
    @test atom_markup(mol.vprops[6]) == [
        [(:default, "C"), (:sub, "2")], [(:default, "H"), (:sub, "7")], [(:default, "N")], [(:sup, "+")]]

    @test false
end

end # model.atom