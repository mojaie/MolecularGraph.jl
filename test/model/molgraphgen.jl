
@testset "model.molgraphgen" begin

@testset "interface" begin
    mol = MolGraphGen()
    @test nv(mol) == 0
    @test ne(mol) == 0
    @test eltype(mol) === Int
    @test edgetype(mol) === Edge{Int}
    @test !is_directed(mol)

    atoms = [SDFAtom(),SDFAtom(),SDFAtom(),SDFAtom(),SDFAtom()]
    bonds = [SDFBond(),SDFBond(),SDFBond(),SDFBond()]
    molstar = MolGraphGen(star_graph(5), atoms, bonds)
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
    nullmol = MolGraphGen()
    @test vproptype(nullmol) === Any
    @test eproptype(nullmol) === Any

    atoms = [SDFAtom(),SDFAtom(),SDFAtom()]
    bonds = [SDFBond(),SDFBond()]
    mol = MolGraphGen(Edge.([(1, 2), (2, 3)]), atoms, bonds, Dict(:hoge => 2))
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
    @test get_prop(mol, 3, 2, :order) == 1
end

@testset "modification" begin
    atoms = [SDFAtom(),SDFAtom(),SDFAtom()]
    bonds = [SDFBond(),SDFBond()]
    mol = MolGraphGen(Edge.([(1, 2), (2, 3)]), atoms, bonds, Dict(:hoge => 2))  # CCC
    @test 3 == set_prop!(mol, :fuga, 3)
    @test get_prop(mol, :fuga) == 3
    @test add_vertex!(mol, SDFAtom(:O))  # CCC.O
    @test get_prop(mol, 4, :symbol) === :O
    @test add_edge!(mol, Edge(3, 4), SDFBond())  # CCCO
    @test degree(mol.graph, 3) == 2
    @test get_prop(mol, Edge(3, 4), :order) == 1
    add_edge!(mol, Edge(1, 4), SDFBond())  # C1CCO1
    @test rem_edge!(mol, Edge(3, 4))  # C1CC.O1
    @test degree(mol.graph, 3) == 1
    @test add_edge!(mol, 4, 3, SDFBond())  # C1CCO1
    @test rem_vertex!(mol, 2)  # COC
    @test degree(mol.graph, 1) == 1
    @test degree(mol.graph, 2) == 2
    @test degree(mol.graph, 3) == 1

    mol = MolGraphGen(smallgraph(:dodecahedral), collect(1:20), collect(1:30))  # MolGraph{Int,Int,Int}
    @test rem_edges!(mol, Edge.([(1, 2), (1, 20), (11, 12), (19, 20)]))
    @test isempty(intersect(mol.eprops, [1, 3, 20, 30]))
    add_edges!(mol, Edge.([(1, 2), (1, 20), (11, 12), (19, 20)]), [1, 3, 20, 30])
    @test [mol.eprops[e] for e in edges(mol)] == collect(1:30)
    vmap = rem_vertices!(mol, [1, 3, 5, 7, 9])
    @test ne(mol) == 16
    @test vmap == [mol.vprops[v] for v in vertices(mol)]
    mol = MolGraphGen(smallgraph(:dodecahedral), collect(1:20), collect(1:30))
    subg, vmap = induced_subgraph(mol, [6, 7, 8, 15, 16])
    @test ne(subg) == 5
    @test vmap == [subg.vprops[v] for v in vertices(subg)]
    subg, vmap = induced_subgraph(mol, Edge.([(1, 2), (1, 11), (1, 20), (4, 20), (19, 20)]))
    @test nv(subg) == 6
    @test sum(values(subg.eprops)) == 45  # edge_rank 1,2,3,9,30
end

@testset "serialization" begin
    atoms = [SDFAtom(),SDFAtom(),SDFAtom()]
    bonds = [SDFBond(),SDFBond()]
    mol = MolGraphGen(Edge.([(1, 2), (2, 3)]), atoms, bonds, Dict(:hoge => 2))
    j = to_json(mol)
    mol2 = MolGraphGen(j)
    @test mol == mol2
    @test mol !== mol2
end

end  # model.molgraph
