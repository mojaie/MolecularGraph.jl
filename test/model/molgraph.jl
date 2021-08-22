
@testset "model.molgraph" begin

@testset "graphmol" begin
    mol = graphmol(SDFileAtom, SDFileBond)
    @test atomcount(mol) == 0
    @test bondcount(mol) == 0

    edges = [(1, 2)]
    atoms = [SDFileAtom(),SDFileAtom(),SDFileAtom()]
    bonds = [SDFileBond()]
    mol = graphmol(edges, atoms, bonds)
    newbond = addbond!(mol, 2, 3, SDFileBond())
    @test newbond == 2
    setatom!(mol, 1, SDFileAtom(:O))
    setbond!(mol, 1, SDFileBond(2))
    @test getatom(mol, 1).symbol == :O
    @test getbond(mol, 1).order == 2

    mol.cache[:test] = [1, 1, 1]
    cp = graphmol(mol)
    @test !haskey(cp.cache, :test)
    cl = clone(mol)
    @test length(cl.cache[:test]) == 3
    @test length(mol.cache[:test]) == 3

    dump = todict(mol)
    loaded = graphmol(dump)
    @test loaded.neighbormap[3][2] == 2
    @test length(loaded.edges) == 2
    @test length(loaded.cache[:test]) == 3
end

end  # model.molgraph
