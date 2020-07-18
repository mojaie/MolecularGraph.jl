
@testset "molgraph" begin

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

@testset "querymol" begin
    q = querymol(SmartsAtom, SmartsBond)
    addatom!(q, SmartsAtom(:Any => true))
    addatom!(q, SmartsAtom(:Any => true))
    addbond!(q, 1, 2, SmartsBond(:Any => true))
    @test atomcount(q) == 2
    @test bondcount(q) == 1
    @test hasbond(q, 1, 2)
end

end  # molgraph
