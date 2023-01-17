
@testset "model.molgraph" begin

@testset "molgraph" begin
    mol = MolGraph()
    @test nv(mol) == 0
    @test ne(mol) == 0
    @test vproptype(mol) === Dict{Symbol,Any}
    @test eproptype(mol) === Dict{Symbol,Any}
    sdfmol = SDFMolGraph()
    @test vproptype(sdfmol) === SDFAtom
    @test eproptype(sdfmol) === SDFBond
    smimol = SMILESMolGraph()
    @test vproptype(smimol) === SMILESAtom
    @test eproptype(smimol) === SMILESBond
    @test !is_directed(smimol)

    @test eltype(SDFMolGraph) === Int
    @test edgetype(SMILESMolGraph) === Edge{Int}
    @test vproptype(SDFMolGraph) === SDFAtom
    @test eproptype(SMILESMolGraph) === SMILESBond

    atoms = [SDFAtom(),SDFAtom(),SDFAtom()]
    bonds = [SDFBond(),SDFBond()]
    moledge = MolGraph(Edge.([(1, 2), (2, 3)]), atoms, bonds)
    @test nv(moledge) == 3
    @test ne(moledge) == 2
    @test vproptype(moledge) === SDFAtom
    @test eproptype(moledge) === SDFBond
    molpg = MolGraph(path_graph(3), atoms, bonds, Dict(:hoge => 2))
    @test nv(molpg) == 3
    @test ne(molpg) == 2
    @test nv(zero(molpg)) == 0
    moliso = MolGraph(Edge.([(1, 2)]), atoms, [SDFBond()])
    @test nv(moliso) == 3
    @test ne(moliso) == 1
end

@testset "access" begin
    atoms = [SDFAtom(),SDFAtom(),SDFAtom(),SDFAtom(),SDFAtom()]
    bonds = [SDFBond(),SDFBond(),SDFBond(),SDFBond()]
    mol = MolGraph(star_graph(5), atoms, bonds)
    @test vertices(mol)[5] == 5
    @test collect(edges(mol))[4] == Edge(1 => 5)
    @test length(neighbors(mol, 1)) == 4
    @test has_vertex(mol, 3)
    @test has_edge(mol, 1, 3)
    @test !has_edge(mol, Edge(2 => 4))
    molgp = MolGraph(Edge.([(1, 2), (2, 3), (2, 4), (3, 5)]), atoms, bonds, Dict(:hoge => 2))
    @test get_prop(molgp, :hoge) == 2
end

@testset "edge_rank" begin
    g = SimpleGraph(Edge.([(3, 2), (6, 5), (1, 2), (3, 4), (5, 4)]))
    @test edge_rank(g, 1, 2) == 1
    @test edge_rank(g, 4, 3) == 3
    @test edge_rank(g, Edge(5, 6)) == 5
end

@testset "edit" begin

end

@testset "serialization" begin

end

end  # model.molgraph
