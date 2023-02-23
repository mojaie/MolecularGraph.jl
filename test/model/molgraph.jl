
@testset "model.molgraph" begin

@testset "interface" begin
    mol = MolGraph()
    @test nv(mol) == 0
    @test ne(mol) == 0
    @test eltype(mol) === Int
    @test edgetype(mol) === Edge{Int}
    @test !is_directed(mol)

    atoms = [SDFAtom(),SDFAtom(),SDFAtom(),SDFAtom(),SDFAtom()]
    bonds = [SDFBond(),SDFBond(),SDFBond(),SDFBond()]
    molstar = MolGraph(star_graph(5), atoms, bonds)
    @test nv(molstar) == 5
    @test ne(molstar) == 4
    @test nv(zero(molstar)) == 0
    @test vertices(molstar)[5] == 5
    @test collect(edges(molstar))[4] == Edge(1 => 5)
    @test length(neighbors(molstar, 1)) == 4
    @test has_vertex(molstar, 3)
    @test has_edge(molstar, 1, 3)
    @test !has_edge(molstar, Edge(2 => 4))
end

@testset "simplemol" begin
    nullmol = MolGraph()
    @test vproptype(nullmol) === Dict{Symbol,Any}
    @test eproptype(nullmol) === Dict{Symbol,Any}
    sdfmol = SDFMolGraph()
    @test vproptype(sdfmol) === SDFAtom
    @test eproptype(sdfmol) === SDFBond
    smimol = SMILESMolGraph()
    @test vproptype(smimol) === SMILESAtom
    @test eproptype(smimol) === SMILESBond

    atoms = [SDFAtom(),SDFAtom(),SDFAtom()]
    bonds = [SDFBond(),SDFBond()]
    mol = MolGraph(Edge.([(1, 2), (2, 3)]), atoms, bonds, Dict(:hoge => 2))
    @test undirectededge(mol, 2, 1) == Edge(1 => 2)
    @test eltype(mol) === Int
    @test edgetype(mol) === Edge{Int}
    @test vproptype(mol) === SDFAtom
    @test eproptype(mol) === SDFBond
    @test length(vprops(mol)) == 3
    @test length(eprops(mol)) == 2
    @test length(props(mol)) == 1
    @test get_prop(mol, :hoge) == 2
    @test !has_prop(mol, :fuga)
    @test get_prop(mol, 1, :symbol) === :C
end

@testset "molgraph" begin
    atoms = [SDFAtom(),SDFAtom(),SDFAtom()]
    bonds = [SDFBond(),SDFBond()]
    mol = MolGraph(Edge.([(1, 2), (2, 3)]), atoms, bonds, Dict(:hoge => 2))
    @test get_prop(mol, 3, 2, :order) == 1
    @test edge_rank(mol, 2, 3) == 2
    @test set_descriptor!(mol, :test, collect(1:3)) == collect(1:3)
    @test has_descriptor(mol, :test)
    @test get_descriptor(mol, :test) == collect(1:3)
end

@testset "editable" begin

end

@testset "serialization" begin

end

end  # model.molgraph
