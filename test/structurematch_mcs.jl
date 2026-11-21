

@testset "structurematch_mcs" begin

@testset "MCS" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    null = smilestomol("")
    @test length(disconnected_mcis(null, null)[1]) == 0
    @test length(disconnected_mces(null, null)[1]) == 0
    @test length(connected_mcis(null, null)[1]) == 0
    @test length(connected_mces(null, null)[1]) == 0
    @test length(tdmcis(null, null)[1]) == 0
    @test length(tdmcis(null, null, tolerance=1)[1]) == 0
    @test length(tdmces(null, null)[1]) == 0
    @test length(tdmces(null, null, tolerance=1)[1]) == 0

    BuOH = smilestomol("CCCO")
    butane = smilestomol("CCCC")
    @test length(disconnected_mcis(butane, BuOH)[1]) == 3
    @test length(disconnected_mces(butane, BuOH)[1]) == 2

    hexane = smilestomol("CCCCCC")
    cyclohexane = smilestomol("C1CCCCC1")
    @test length(tdmcis(hexane, cyclohexane)[1]) == 4
    @test length(tdmcis(hexane, cyclohexane, tolerance=1)[1]) == 4
    @test length(tdmces(hexane, cyclohexane)[1]) == 4
    @test length(tdmces(hexane, cyclohexane, tolerance=1)[1]) == 4

    tms = smilestomol("C[Si](C)(C)C")
    tsm = smilestomol("[SiH3]C([SiH3])([SiH3])[SiH3]")  # bracket atoms have no implicit hydrogens
    @test length(disconnected_mcis(tms, tsm)[1]) == 2
    @test length(disconnected_mces(tms, tsm)[1]) == 1
    @test length(tdmcis(tms, tsm)[1]) == 2
    @test length(tdmces(tms, tsm)[1]) == 1

    sulfide = smilestomol("CSSC")
    disconn = smilestomol("CS.SC")
    @test length(connected_mcis(sulfide, disconn)[1]) == 2
    @test length(disconnected_mcis(sulfide, disconn)[1]) == 3
    @test length(connected_mces(sulfide, disconn)[1]) == 1
    @test length(disconnected_mces(sulfide, disconn)[1]) == 2
    @test length(tdmcis(sulfide, disconn, tolerance=1)[1]) == 2
    @test length(tdmces(sulfide, disconn, tolerance=1)[1]) == 1

    # Delta-Y
    cyclopropane = smilestomol("C1CC1")
    isopropane = smilestomol("CC(C)C")
    @test length(disconnected_mcis(cyclopropane, isopropane)[1]) == 2
    @test length(disconnected_mces(cyclopropane, isopropane)[1]) == 2

    tetrahedrane = smilestomol("C12C3C1C23")
    fused = smilestomol("C1C2C1C2")
    spiro = smilestomol("C1CC12CC2")
    @test length(disconnected_mcis(tetrahedrane, fused)[1]) == 3
    @test length(disconnected_mcis(tetrahedrane, spiro)[1]) == 3
    @test length(disconnected_mces(tetrahedrane, fused)[1]) == 5
    @test length(disconnected_mces(tetrahedrane, spiro)[1]) == 4

    cid6437877 = smilestomol("CC1=C(SC=N1)C=CC2=C(N3C(C(C3=O)NC(=O)C(=NOC)C4=CSC(=N4)N)SC2)C(=O)OCOC(=O)C(C)(C)C")  # Cefditoren pivoxil
    cid5481173 = smilestomol("CC(C)(C(=O)O)ON=C(C1=CSC(=N1)N)C(=O)NC2C3N(C2=O)C(=C(CS3)C[N+]4=CC=CC=C4)C(=O)[O-]")  # Ceftazidime
    @test length(tdmcis(cid6437877, cid5481173, tolerance=1)[1]) == 28
    @test length(tdmces(cid6437877, cid5481173, tolerance=1)[1]) == 29
    @test length(tdmcis(cid6437877, cid5481173, diameter=8)[1]) == 18
    @test length(tdmces(cid6437877, cid5481173, diameter=8)[1]) == 20
    # global_logger(default_logger)
end

@testset "mcs_constraints" begin
    BuOH = smilestomol("CCCO")
    vvec = vmatchvecgen(BuOH)
    @test [vvec(i) for i in 1:4] == [366, 366, 366, 368]  # sp3 * 120 + atom number
    @test ematchvecgen(BuOH)(1, 2) == 0

    c = mcis_constraints(BuOH)
    @test c isa ConstraintArrayMCIS{Int,Bool,Int,Int}
    @test c.nv == 4
    @test c.pairs == [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    @test c.distances == Bool[1, 0, 0, 1, 0, 1]
    @test c.vattrs1 == [vvec(u) for (u, v) in c.pairs]
    @test c.vattrs2 == [vvec(v) for (u, v) in c.pairs]
    c = mces_constraints(BuOH)
    @test c isa ConstraintArrayMCES{Int,Bool,Int,Int}
    @test c.nv == 3  # bonds
    @test issetequal(values(c.revmap), edges(BuOH.graph))

    hexane = smilestomol("CCCCCC")
    c = tdmcis_constraints(hexane)
    @test c isa ConstraintArrayMCIS{Int,UInt8,Int,Int}
    @test length(c.pairs) == 15
    @test maximum(c.distances) == 5
    c = tdmcis_constraints(hexane, diameter=2)
    @test length(c.pairs) == 9
    @test maximum(c.distances) == 2
    c = tdmces_constraints(hexane)
    @test c isa ConstraintArrayMCES{Int,UInt8,Int,Int}
    @test c.nv == 5

    # serialization
    cyclohexane = smilestomol("C1CCCCC1")
    for constraints in (tdmcis_constraints, tdmces_constraints)
        h = constraints(hexane)
        hj = JSON.parse(JSON.json(h), typeof(h))
        @test hj.pairs == h.pairs
        @test hj.distances == h.distances
        @test length(maximum_common_subgraph(constraints(cyclohexane), hj)[1]) == 4
    end
end

@testset "mcs_custom_match" begin
    BuOH = smilestomol("CCCO")
    butane = smilestomol("CCCC")
    anyatom = mol -> (i -> 0)
    @test length(disconnected_mcis(butane, BuOH)[1]) == 3
    @test length(maximum_common_subgraph(
        mcis_constraints(butane, anyatom), mcis_constraints(BuOH, anyatom))[1]) == 4

    function bondorder(mol)
        ernk = edge_rank(mol)
        bo = bond_order(mol)
        return (u, v) -> has_edge(mol.graph, u, v) ? bo[edge_rank(ernk, u, v)] : 0
    end
    propene = smilestomol("C=CC")
    propane = smilestomol("CCC")
    @test length(maximum_common_subgraph(
        mces_constraints(propene, anyatom), mces_constraints(propane, anyatom))[1]) == 2
    @test length(maximum_common_subgraph(
        mces_constraints(propene, anyatom, bondorder),
        mces_constraints(propane, anyatom, bondorder))[1]) == 1
    @test length(maximum_common_subgraph(
        mcis_constraints(propene, anyatom, bondorder),
        mcis_constraints(propane, anyatom, bondorder))[1]) == 2
    @test length(tdmcis(propene, propane, vmatchvecgen=anyatom)[1]) == 3
    @test length(tdmcis(propene, propane, vmatchvecgen=anyatom, ematchvecgen=bondorder)[1]) == 2
end

@testset "mcs_mapping" begin
    EtOH = smilestomol("CCO")
    @test disconnected_mcis(EtOH, smilestomol("OCC")) == (Dict(1 => 3, 2 => 2, 3 => 1), :done)
    emap, status = disconnected_mces(EtOH, smilestomol("OCC"))
    @test status == :done
    @test emap == Dict(Edge(1, 2) => Edge(2, 3), Edge(2, 3) => Edge(1, 2))

    # mapped atoms have the same symbol and the connected MCS is connected
    mol1 = smilestomol("CC1=C(SC=N1)C=CC2=C(N3C(C(C3=O)NC(=O)C(=NOC)C4=CSC(=N4)N)SC2)C(=O)OCOC(=O)C(C)(C)C")
    mol2 = smilestomol("CC(C)(C(=O)O)ON=C(C1=CSC(=N1)N)C(=O)NC2C3N(C2=O)C(=C(CS3)C[N+]4=CC=CC=C4)C(=O)[O-]")
    sym1 = atom_symbol(mol1)
    sym2 = atom_symbol(mol2)
    nmap, status = connected_mcis(mol1, mol2)
    @test status == :done
    @test length(nmap) == 25
    @test all(sym1[u] == sym2[v] for (u, v) in nmap)
    @test is_connected(induced_subgraph(mol1.graph, collect(keys(nmap)))[1])
    @test is_connected(induced_subgraph(mol2.graph, collect(values(nmap)))[1])
    emap, status = connected_mces(mol1, mol2)
    @test status == :done
    @test length(emap) == 27
    @test issubset(keys(emap), edges(mol1.graph))
    @test issubset(values(emap), edges(mol2.graph))
    @test is_connected(induced_subgraph(mol1.graph, collect(keys(emap)))[1])
end

@testset "mcs_timeout_targetsize" begin
    # Cefditoren pivoxil and Ceftazidime: disconnected MCS takes much longer than 1 sec
    mol1 = smilestomol("CC1=C(SC=N1)C=CC2=C(N3C(C(C3=O)NC(=O)C(=NOC)C4=CSC(=N4)N)SC2)C(=O)OCOC(=O)C(C)(C)C")
    mol2 = smilestomol("CC(C)(C(=O)O)ON=C(C1=CSC(=N1)N)C(=O)NC2C3N(C2=O)C(=C(CS3)C[N+]4=CC=CC=C4)C(=O)[O-]")
    for f in (disconnected_mcis, disconnected_mces, connected_mcis, connected_mces, tdmcis, tdmces)
        res, status = f(mol1, mol2, timeout=0)
        @test status == :timedout
        res, status = f(mol1, mol2, targetsize=8)
        @test status == :targetreached
        @test length(res) == 8
    end
    # suboptimal results within the time limit
    for f in (disconnected_mcis, disconnected_mces)
        elapsed = @elapsed res, status = f(mol1, mol2, timeout=0.5)
        @test status == :timedout
        @test !isempty(res)
        @test elapsed < 5  # default timeout is 10 sec
    end
    # tolerance is passed together with the options
    res, status = tdmcis(mol1, mol2, tolerance=1, targetsize=10)
    @test status == :targetreached
    @test length(res) == 10
end

end # structurematch_mcs
