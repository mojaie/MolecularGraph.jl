
@testset "querycontainment" begin

@testset "truthtable" begin
    # truthtable generation
    anytrue = QueryAtom(Tuple{Int,Int}[], [qanytrue()])
    tt1, tt2 = generate_truthtable(anytrue, anytrue)
    @test isempty(tt1.props)
    @test isempty(tt2.props)
    @test tt1.func([])
    @test tt2.func([])

    # query match
    op1 = QueryAtom(
        [(1, 2), (1, 3), (1, 5), (3, 4), (5, 6), (5, 7)],
        [qand(), qeq(:symbol, "C"), qnot(), qtrue(:isaromatic), qor(),
        qeq(:connectivity, "3"), qeq(:connectivity, "4")])
    propmap = querypropmap(op1)
    @test issetequal(keys(propmap), [:symbol, :isaromatic, :connectivity])
    tbl1 = QueryTruthTable(op1)
    func1 = QueryTruthTable(
        v -> v[4] & ~v[3] & (v[1] | v[2]),
        [qeq(:connectivity, "3"), qeq(:connectivity, "4"),
        qtrue(:isaromatic), qeq(:symbol, "C")]
    )
    @test querymatch(tbl1, func1, true)

    # tautology
    """
    q \\ r  T other F
    T      T   F   F
    other  T   ?   F
    F      T   T   T
    """
    props = [qeq(:symbol, "C")]
    t = QueryTruthTable(v -> true, props)
    f = QueryTruthTable(v -> false, props)
    o = QueryTruthTable(v -> v[1], props)
    @test issubset(t, t)
    @test issubset(f, f)
    @test !issubset(t, o)
    @test !issubset(t, f)
    @test issubset(o, t)
    @test !issubset(o, f)
    @test issubset(f, t)
    @test issubset(f, o)

    # containment
    props = [qeq(:symbol, "C"), qeq(:symbol, "N")]
    symC = QueryTruthTable(v -> v[1], props)
    symN = QueryTruthTable(v -> v[2], props)
    notN = QueryTruthTable(v -> ~v[2], props)
    props = [qeq(:total_hydrogens, "1"), qeq(:charge, "1")]
    hc1 = QueryTruthTable(v -> v[1], props)
    ch1 = QueryTruthTable(v -> v[2], props)
    @test symC == symC
    @test symC != symN
    @test symN != symC
    @test symN != notN
    @test notN != symN
    @test hc1 == hc1
    @test hc1 != ch1
    @test ch1 != hc1

    props = [qeq(:charge, "1"), qtrue(:isaromatic), qeq(:symbol, "N")]
    and1 = QueryTruthTable(v -> v[2] & v[3], props)
    and2 = QueryTruthTable(v -> v[1] & v[2] & v[3], props)
    @test issubset(and2, and1)
    @test !issubset(and1, and2)

    props = [qeq(:charge, "-1"), qeq(:charge, "0"), qeq(:charge, "1")]
    or1 = QueryTruthTable(v -> v[2] | v[3], props)
    or2 = QueryTruthTable(v -> v[1] | v[2] | v[3], props)
    @test issubset(or1, or2)
    @test !issubset(or2, or1)

    props = [
        qeq(:charge, "0"), qeq(:charge, "1"), qtrue(:isaromatic),
        qeq(:smallest_ring, "6"), qeq(:symbol, "N"), qeq(:symbol, "O"), qeq(:symbol, "S")
    ]
    nested1 = QueryTruthTable(
        v -> ~v[3] & (v[5] | v[6]) & (v[1] | v[2]), props)
    nested2 = QueryTruthTable(
        v -> ~v[3] & (v[5] | v[6] | v[7]) & (v[1] | v[2] | v[4]), props)
    nested3 = QueryTruthTable(
        v -> ~v[3] & v[6] & (v[1] | v[2]), props)
    @test issubset(nested1, nested2)
    @test !issubset(nested2, nested1)
    @test issubset(nested3, nested1)
    @test !issubset(nested1, nested3)
end

@testset "resolve_disjoint_not" begin
    oors = QueryAtom(
        [(1, 2), (1, 3)],
        [qor(), qeq(:symbol, "O"), qeq(:symbol, "S")])  # [#8,#16]
    notn = QueryAtom([(1, 2)], [qnot(), qeq(:symbol, "N")])  # [!#7]
    oors_ = deepcopy(oors)
    resolve_disjoint_not!(oors_, querypropmap(notn))  # [#8,#16]
    @test oors == oors_
    resolve_disjoint_not!(notn, querypropmap(oors))  # [!#7,#8,#16]
    @test notn == QueryAtom(
        [(1, 2), (1, 3), (1, 4), (4, 5)],
        [qor(), qeq(:symbol, "O"), qeq(:symbol, "S"),
        qnot(), qeq(:symbol, "N")])
    noto = QueryAtom([(1, 2)], [qnot(), qeq(:symbol, "O")])  # [!#8]
    resolve_disjoint_not!(noto, querypropmap(oors))  # [!#8,#16]
    @test noto == QueryAtom(
        [(1, 2), (1, 3), (3, 4)],
        [qor(), qeq(:symbol, "S"), qnot(), qeq(:symbol, "O")])
    notarom = QueryAtom([(1, 2)], [qnot(), qtrue(:isaromatic)])  # [a]
    notarom_ = deepcopy(notarom)
    resolve_disjoint_not!(notarom, querypropmap(oors))  # [a]
    @test notarom_ == notarom
end

@testset "resolve_recursive" begin
    rec1 = QueryAtom(Tuple{Int,Int}[], [qeq(:recursive, "[#6][NH]")])  # [$([#6][NH])]
    rec2 = QueryAtom(Tuple{Int,Int}[], [qeq(:recursive, "[#6]N")])  # [$([#6]N)]
    rec1_ = deepcopy(rec1)
    resolve_recursive!(rec1_, querypropmap(rec2))  # [$([#6][NH]);#6;$([#6]N)]
    @test rec1_ == QueryAtom(
        [(1, 2), (1, 3), (1, 4)],
        [qand(), qeq(:symbol, "C"), qeq(:recursive, "[#6][NH]"), qeq(:recursive, "[#6]N")])
    rec2_ = deepcopy(rec2)
    resolve_recursive!(rec2_, querypropmap(rec1))  # [$([#6]N);#6]
    @test rec2_ == QueryAtom(
        [(1, 2), (1, 3)],
        [qand(), qeq(:symbol, "C"), qeq(:recursive, "[#6]N")])

    nc = QueryAtom(Tuple{Int,Int}[], [qeq(:recursive, "[nH][#6]")]) # [$([nH][#6])]
    aromn = QueryAtom(
        [(1, 2), (1, 3), (1, 4)],
        [qand(), qeq(:symbol, "N"), qtrue(:isaromatic), qeq(:total_hydrogens, "1")]) # [nH]
    nonan = QueryAtom(
        [(1, 2), (1, 3), (3, 4)],
        [qand(), qeq(:symbol, "N"), qnot(), qtrue(:isaromatic)]) # N
    nc_ = deepcopy(nc)
    resolve_recursive!(nc_, querypropmap(aromn))  # [$([nH][#6]);[nH]]
    @test nc_ == QueryAtom(
        [(1, 2), (1, 3), (3, 4), (3, 5), (5, 6), (5, 7)],
        [qand(), qeq(:recursive, "[nH][#6]"),
        qand(), qeq(:total_hydrogens, "1"), qand(), qeq(:symbol, "N"), qtrue(:isaromatic)])
    nc2_ = deepcopy(nc)
    resolve_recursive!(nc2_, querypropmap(nonan))  # [$([nH][#6]);[nH]]
    @test nc_ == nc2_
    aromn_ = deepcopy(aromn)
    resolve_recursive!(aromn_, querypropmap(nonan))  # [nH]
    @test aromn_ == aromn
    nonan_ = deepcopy(nonan)
    resolve_recursive!(nonan_, querypropmap(aromn))  # N
    @test nonan_ == nonan

    diazo1 = QueryAtom(Tuple{Int,Int}[], [qeq(:recursive, "N=N=C")])
    diazo2 = QueryAtom(
        [(1, 2), (1, 3)],
        [qor(), qeq(:recursive, "N=N=C"), qeq(:recursive, "N#N-C")])
    diazo1_ = deepcopy(diazo1)
    resolve_recursive!(diazo1_, querypropmap(diazo2))
    @test diazo1_ == QueryAtom(
        [(1, 2), (2, 3), (2, 4), (4, 5), (1, 6)],
        [qand(), qand(), qeq(:symbol, "N"), qnot(), qtrue(:isaromatic),
        qeq(:recursive, "N=N=C")])
    diazo2_ = deepcopy(diazo2)
    resolve_recursive!(diazo2_, querypropmap(diazo1))
    @test diazo2_ == QueryAtom(
        [(1, 2), (1, 3), (2, 4), (2, 5), (5, 6), (5, 7), (7, 8),
        (3, 9), (3, 10), (10, 11), (10, 12), (12, 13)],
        [qor(), qand(), qand(),
        qeq(:recursive, "N=N=C"), qand(), qeq(:symbol, "N"), qnot(), qtrue(:isaromatic),
        qeq(:recursive, "N#N-C"), qand(), qeq(:symbol, "N"), qnot(), qtrue(:isaromatic)])
end

@testset "smarts" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))

    # no properties
    fmr = smartstomol("*1***1")
    three = smartstomol("***")
    @test has_substruct_match(fmr, three)
    @test !has_substruct_match(three, fmr)

    tmse = smartstomol("O[Si](C)(C)C")
    tms = smartstomol("C[Si](C)C")
    @test has_substruct_match(tmse, tms)
    @test !has_exact_match(tmse, tms)

    primary = smartstomol("[CH2][OH]")
    alcohol = smartstomol("C[OH]")
    @test has_exact_match(primary, alcohol)
    @test !has_exact_match(alcohol, primary)

    pyrrole = smartstomol("[nH]1cccc1")
    pyridine = smartstomol("n1ccccc1")
    pyrrolidine = smartstomol("N1CCCC1")
    cyclicamine = smartstomol("[#7]1[#6][#6][#6][#6]1")
    @test !has_substruct_match(pyridine, pyrrole)
    @test !has_substruct_match(pyrrole, pyrrolidine)
    @test has_exact_match(pyrrole, cyclicamine)
    @test has_exact_match(pyrrolidine, cyclicamine)

    carbonylazide = smartstomol("O=CN=[N+]=[N-]")
    azide = smartstomol("N=[N+]=[N-]")
    hydrazine = smartstomol("N=N")
    stricthydrazine = smartstomol("[N+0]=[N+0]")
    @test has_substruct_match(carbonylazide, azide)
    @test has_substruct_match(azide, hydrazine)
    @test !has_substruct_match(azide, stricthydrazine)

    halo = smartstomol("[#9,#17,#35]c1ccccc1")
    narrow = smartstomol("[#9,#17]c1ccccc1")
    broad = smartstomol("[#9,#17,#35,#53]c1ccccc1")
    broad2 = smartstomol("[#9,35#17,#35,#53]c1ccccc1")
    @test has_exact_match(narrow, halo)
    @test !has_exact_match(broad, halo)
    @test has_exact_match(halo, broad)
    @test !has_exact_match(halo, broad2)

    hetero = smartstomol("[!#6&!#7&!#8]1[#6][#6][#6][#6][#6]1")
    oxo = smartstomol("[#8]1[#6][#6][#6][#6][#6]1")
    sulfo = smartstomol("[#16]1[#6][#6][#6][#6][#6]1")
    thiopyrylium = smartstomol("[s+]1ccccc1")
    @test !has_exact_match(oxo, hetero)
    @test has_exact_match(sulfo, hetero)
    @test has_exact_match(thiopyrylium, sulfo)
    @test has_exact_match(thiopyrylium, hetero)

    diazocarbonyl1 = smartstomol(raw"[$(N=N=C~C=O)]")
    diazocarbonyl2 = smartstomol(raw"[$(N=N=C~C=O),$(N#N-C~C=O)]")
    @test has_exact_match(diazocarbonyl1, diazocarbonyl2)
    @test !has_exact_match(diazocarbonyl2, diazocarbonyl1)

    nested = smartstomol(raw"[$([CH]=[$(NOC)])]C=O")
    nestedor = smartstomol(raw"[$([CH]=[$(NOC),$(NO[Si])])]")
    @test has_substruct_match(nested, nestedor)
    @test !has_substruct_match(nestedor, nested)

    naryl = smartstomol(raw"OP(=O)(=[S,O])[$(Na)]")
    phos = smartstomol(raw"OP(=O)(=[S,O])N")
    phos2 = smartstomol(raw"[$([S,O]=PN)]")
    @test has_exact_match(naryl, phos)
    @test !has_exact_match(phos, naryl)
    @test !has_substruct_match(phos, phos2)

    tfas = smartstomol(raw"C(F)(F)(F)C(=O)S")
    tfmk1 = smartstomol(raw"PS[$(C(=O))](=O)C(F)(F)(F)")
    tfmk2 = smartstomol(raw"[$(C(=O));!$(C-N);!$(C-O);!$(C-S)]C(F)(F)(F)")
    @test has_substruct_match(tfmk1, tfas)
    @test !has_substruct_match(tfas, tfmk2)

    quart = smartstomol(raw"[C+,Cl+,I+,P+,S+]")
    sulfonium = smartstomol(raw"[S+;X3;$(S-C);!$(S-[O;D1])]")
    cys = smartstomol(raw"NC(C=O)CS")
    @test has_exact_match(sulfonium, quart)
    @test !has_substruct_match(cys, sulfonium)
    @test !has_substruct_match(cys, quart)

    halobenzene = smartstomol(raw"c1c([O;D1])c(-[Cl,Br,I])cc(-[Cl,Br,I])c1")
    halo = smartstomol(raw"c[F,Cl,Br,I]")
    @test has_substruct_match(halobenzene, halo)

    pivoxil = smartstomol(raw"OCOC(=O)C([CH3])([CH3])[CH3]")
    ester = smartstomol(raw"[#6]C(=O)O[#6]")
    ester2 = smartstomol(raw"[#6]-C(=O)O-[#6]")  # MLSMR structural alert Ester
    @test has_substruct_match(pivoxil, ester)
    @test !has_substruct_match(pivoxil, ester2)

    carbon = smartstomol("C")
    # BMS structural alert 'contains metal'
    # Just for int overflow check (set max property size cutoff, do not generate 2^61 truthtable)
    metals = smartstomol(raw"[$([Ru]),$([Rh]),$([Se]),$([Pd]),$([Sc]),$([Bi]),$([Sb]),$([Ag]),$([Ti]),$([Al]),$([Cd]),$([V]),$([In]),$([Cr]),$([Sn]),$([Mn]),$([La]),$([Fe]),$([Er]),$([Tm]),$([Yb]),$([Lu]),$([Hf]),$([Ta]),$([W]),$([Re]),$([Co]),$([Os]),$([Ni]),$([Ir]),$([Cu]),$([Zn]),$([Ga]),$([Ge]),$([As]),$([Y]),$([Zr]),$([Nb]),$([Ce]),$([Pr]),$([Nd]),$([Sm]),$([Eu]),$([Gd]),$([Tb]),$([Dy]),$([Ho]),$([Pt]),$([Au]),$([Hg]),$([Tl]),$([Pb]),$([Ac]),$([Th]),$([Pa]),$([Mo]),$([U]),$([Tc]),$([Te]),$([Po]),$([At])]")
    @test !has_substruct_match(carbon, metals)

    # global_logger(default_logger)
end

end