
@testset "model.molgraph" begin

@testset "interface" begin
    mol = MolGraph()
    @test nv(mol) == 0
    @test ne(mol) == 0
    @test eltype(mol) === Int
    @test edgetype(mol) === Edge{Int}
    @test !is_directed(mol)

    fe = MolGraph(Edge{Int}[], [SDFAtom(:Fe)], SDFBond[])
    @test nv(fe) == 1
    @test ne(fe) == 0

    atoms = [SDFAtom(),SDFAtom(),SDFAtom(),SDFAtom(),SDFAtom()]
    bonds = [SDFBond(),SDFBond(),SDFBond(),SDFBond()]
    molstar = MolGraph(collect(edges(star_graph(5))), atoms, bonds)
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
    @test vproptype(nullmol) === Any
    @test eproptype(nullmol) === Any
    sdfmol = SDFMolGraph()
    @test vproptype(sdfmol) === SDFAtom
    @test eproptype(sdfmol) === SDFBond
    smimol = SMILESMolGraph()
    @test vproptype(smimol) === SMILESAtom
    @test eproptype(smimol) === SMILESBond

    atoms = [SDFAtom(),SDFAtom(),SDFAtom()]
    bonds = [SDFBond(),SDFBond()]
    mol = MolGraph(Edge.([(1, 2), (2, 3)]), atoms, bonds, gprop_map=Dict(:hoge => 2))
    @test u_edge(mol, 2, 1) == Edge(1 => 2)
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
    mol = MolGraph(Edge.([(1, 2), (2, 3)]), atoms, bonds, gprop_map=Dict(:hoge => 2))
    @test get_prop(mol, 3, 2, :order) == 1
    @test edge_rank(mol, 2, 3) == 2
    @test set_state!(mol, :test, collect(1:3)) == collect(1:3)
    @test has_state(mol, :test)
    @test get_state(mol, :test) == collect(1:3)
end

@testset "modification" begin
    atoms = [SDFAtom(),SDFAtom(),SDFAtom()]
    bonds = [SDFBond(),SDFBond()]
    mol = MolGraph(Edge.([(1, 2), (2, 3)]), atoms, bonds, gprop_map=Dict(:hoge => 2))  # CCC
    set_prop!(mol, :fuga, 3)
    @test get_prop(mol, :fuga) == 3
    @test add_vertex!(mol, SDFAtom(:O))  # CCC.O
    @test get_prop(mol, 4, :symbol) === :O
    @test add_edge!(mol, Edge(3, 4), SDFBond())  # CCCO
    @test degree(mol.graph, 3) == 2
    @test get_prop(mol, Edge(3, 4), :order) == 1
    # @test length(mol.edge_rank) == 3    move to properties
    add_edge!(mol, Edge(1, 4), SDFBond())  # C1CCO1
    @test rem_edge!(mol, Edge(3, 4))  # C1CC.O1
    @test degree(mol.graph, 3) == 1
    # @test length(mol.edge_rank) == 3    move to properties
    @test add_edge!(mol, 4, 3, SDFBond())  # C1CCO1
    @test rem_vertex!(mol, 2)  # COC
    @test degree(mol.graph, 1) == 1
    @test degree(mol.graph, 2) == 2
    @test degree(mol.graph, 3) == 1
    @test length(mol.edge_rank) == 2

    mol = MolGraph(collect(edges(smallgraph(:dodecahedral))), collect(1:20), collect(1:30))  # MolGraph{Int,Int,Int}
    for e in Edge.([(1, 2), (1, 20), (11, 12), (19, 20)])
        rem_edge!(mol, e)
    end
    @test isempty(intersect(mol.eprops, [1, 3, 20, 30]))
    for (e, p) in zip(Edge.([(1, 2), (1, 20), (11, 12), (19, 20)]), [1, 3, 20, 30])
        add_edge!(mol, e, p)
    end
    @test issetequal(values(mol.eprops), collect(1:30))
    vmap = rem_vertices!(mol, [1, 3, 5, 7, 9])
    @test ne(mol) == 16
    @test issetequal(vmap, values(mol.vprops))
    mol = MolGraph(collect(edges(smallgraph(:dodecahedral))), collect(1:20), collect(1:30))
    subg, vmap = induced_subgraph(mol, [6, 7, 8, 15, 16])
    @test ne(subg) == 5
    @test issetequal(vmap, values(subg.vprops))
    subg, vmap = induced_subgraph(mol, Edge.([(1, 2), (1, 11), (1, 20), (4, 20), (19, 20)]))
    @test nv(subg) == 6
    @test sum(values(subg.eprops)) == 45  # edge_rank 1,2,3,9,30
    rem_vertex!(mol, 20)
    @test nv(mol) == 19
    rem_vertices!(mol, [17, 18, 19])
    @test nv(mol) == 16
end

@testset "serialization" begin
    atoms = [SDFAtom(),SDFAtom(),SDFAtom()]
    bonds = [SDFBond(),SDFBond()]
    mol = MolGraph(Edge.([(1, 2), (2, 3)]), atoms, bonds, gprop_map=Dict(:hoge => 2))
    mol2 = MolGraph(to_json(mol))
    @test mol == mol2
    @test mol !== mol2

    atoms = [SMILESAtom(),SMILESAtom(),SMILESAtom()]
    bonds = [SMILESBond(),SMILESBond()]
    mol = MolGraph(Edge.([(1, 2), (2, 3)]), atoms, bonds, gprop_map=Dict(:hoge => 2))
    mol2 = MolGraph(to_json(mol))
    @test mol == mol2
    @test mol !== mol2

    mol = MolGraph()
    mol2 = MolGraph(to_json(mol))
    @test mol == mol2
    @test mol !== mol2

    atoms = [SDFAtom(),SDFAtom(),SDFAtom()]
    bonds = SDFBond[]
    mol = MolGraph(Edge{Int}[], atoms, bonds)
    mol2 = MolGraph(to_json(mol))
    @test mol == mol2
    @test mol !== mol2
end

end  # model.molgraph
