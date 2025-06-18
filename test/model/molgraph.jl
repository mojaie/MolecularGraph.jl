
@testset "model.molgraph" begin

@testset "graph_interface" begin
    nullmol = MolGraph()
    nullsdf = SDFMolGraph()
    nullsmiles = SMILESMolGraph()
    fe = MolGraph(Edge{Int}[], [SDFAtom(;symbol=:Fe)], SDFBond[])
    atoms = [SDFAtom(), [SDFAtom(symbol=:H) for _ in 1:4]...]
    bonds = [SDFBond() for i in 1:4]
    methane = MolGraph(collect(edges(star_graph(5))), atoms, bonds)
    atoms = [SDFAtom(), SDFAtom(), [SDFAtom(symbol=:H) for _ in 1:6]...]
    bonds = [SDFBond() for i in 1:7]
    ethane = MolGraph(
        Edge.([(1, 2), (1, 3), (1, 4), (1, 5), (2, 6), (2, 7), (2, 8)]),
        atoms, bonds)

    # compatibility with Graphs.jl
    @test eltype(nullmol) === Int
    @test edgetype(nullmol) === Edge{Int}
    @test !is_directed(nullmol)
    @test nv(nullmol) == 0
    @test nv(fe) == 1
    @test nv(methane) == 5
    @test nv(zero(methane)) == 0
    @test ne(nullmol) == 0
    @test ne(fe) == 0
    @test ne(methane) == 4
    @test vertices(methane)[5] == 5
    @test collect(edges(methane))[4] == Edge(1 => 5)
    @test length(neighbors(methane, 1)) == 4
    @test has_vertex(methane, 3)
    @test has_edge(methane, 1, 3)
    @test !has_edge(methane, Edge(2 => 4))

    # custom graph methods
    @test edge_rank(ethane, 2, 7) == 6
    @test u_edge(methane, 2, 1) == u_edge(methane, 1, 2)
    @test u_edge(methane, 4, 1) == u_edge(methane, Edge(1 => 4))
    @test u_edge(methane, 3, 1) == u_edge(eltype(edgetype(methane)), 1, 3)
    @test length(ordered_neighbors(methane, 1)) == 4
    @test length(edge_neighbors(ethane, 1, 2)[1]) == 3
    @test length(ordered_edge_neighbors(ethane, Edge(1 => 2))[1]) == 3

    # node, edge or graph properties
    @test vproptype(nullmol) === SDFAtom
    @test vproptype(nullsdf) === SDFAtom
    @test vproptype(nullsmiles) === SMILESAtom
    @test eproptype(nullmol) === SDFAtom
    @test eproptype(nullsdf) === SDFBond
    @test eproptype(nullsmiles) === SMILESBond
    @test get_prop(methane, 1, :symbol) === :C
    @test get_prop(methane, 1, 2, :order) === 1
end


@testset "modification" begin
    atoms = [SDFAtom(),SDFAtom(),SDFAtom()]
    bonds = [SDFBond(),SDFBond()]
    mol = MolGraph(Edge.([(1, 2), (2, 3)]), atoms, bonds)  # CCC
    @test !mol.state.has_updates  # dispatch_update! called after initialization
    # edit graph properties
    @test add_vertex!(mol, SDFAtom(;symbol=:O))  # CCC.O
    @test mol.state.has_updates
    mol.state.has_updates = false
    @test get_prop(mol, 4, :symbol) === :O

    @test add_edge!(mol, Edge(3, 4), SDFBond())  # CCCO
    @test mol.state.has_updates
    mol.state.has_updates = false
    @test degree(mol.graph, 3) == 2
    @test get_prop(mol, Edge(3, 4), :order) == 1

    add_edge!(mol, Edge(1, 4), SDFBond())  # C1CCO1
    @test mol.state.has_updates
    mol.state.has_updates = false

    @test rem_edge!(mol, Edge(3, 4))  # C1CC.O1
    @test mol.state.has_updates
    mol.state.has_updates = false
    @test degree(mol.graph, 3) == 1
    @test issetequal(keys(mol.eprops), Edge.([(1, 2), (1, 4), (2, 3)]))

    @test add_edge!(mol, 4, 3, SDFBond())  # C1CCO1
    @test rem_vertex!(mol, 2)  # COC
    @test mol.state.has_updates
    mol.state.has_updates = false
    @test degree(mol.graph, 1) == 1
    @test degree(mol.graph, 2) == 2
    @test degree(mol.graph, 3) == 1
    @test issetequal(keys(mol.vprops), collect(1:3))
    @test issetequal(keys(mol.eprops), Edge.([(1, 2), (2, 3)]))
    @test rem_vertex!(mol, 2)  # C.C

    # Vertex re-ordering when vertices/edges are removed
    mol = MolGraph(
        collect(edges(smallgraph(:dodecahedral))), collect(1:20), collect(1:30)
    )  # MolGraph{Int,Int,Int}
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
    vmap = rem_vertices!(mol, collect(1:14))
    @test ne(mol) == 0
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


@testset "metadata" begin
    gprop = MolProperty{Int}()
    data = Dict(
        "id1" => "hoge", "id2" => "fuga", "valid" => "true", "exp_pka" => "9.24", "stock_mg" => "25"
    )
    gprop.metadata = reconstruct(Val(:metadata), MolProperty{Int}, data)
    atoms = [SDFAtom(),SDFAtom(),SDFAtom()]
    bonds = [SDFBond(),SDFBond()]
    mol = MolGraph(Edge.([(1, 2), (2, 3)]), atoms, bonds, gprops=gprop)
    # getter/setter
    @test get_prop(mol, "exp_pka") == "9.24"
    set_prop!(mol, "new_id", "CP00001")
    @test get_prop(mol, "new_id") == "CP00001"
    @test !has_prop(mol, "unknown")
    # comvenient methods
    @test mol["valid"] == "true"
    mol["deleterious"] = "false"
    @test mol["deleterious"] == "false"
end

end  # model.molgraph
