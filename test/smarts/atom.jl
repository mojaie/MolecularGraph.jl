
@testset "smarts.atom" begin

@testset "atomsymbol" begin
    # outside of []
    state = SMARTSParser{SMARTSMolGraph}("Cl")
    qtree = QueryAtom()
    cl = atomsymbol!(state, qtree)
    @test cl == 1
    @test qtree == QueryAtom(Tuple{Int,Int}[], [qeq(:symbol, "Cl")])

    state = SMARTSParser{SMARTSMolGraph}("Cr")  # only Cr inside [] can be recognized
    qtree = QueryAtom()
    cr = atomsymbol!(state, qtree)
    @test cr == 1
    @test qtree == QueryAtom(
        [(1, 2), (1, 3), (3, 4)], [qand(), qeq(:symbol, "C"), qnot(), qtrue(:isaromatic)])
    @test state.pos == 2  # Next, read 'r'

    state = SMARTSParser{SMARTSMolGraph}("p")
    qtree = QueryAtom()
    aromp = atomsymbol!(state, qtree)
    @test aromp == 1
    @test qtree == QueryAtom(
        [(1, 2), (1, 3)], [qand(), qeq(:symbol, "P"), qtrue(:isaromatic)])

    state = SMARTSParser{SMARTSMolGraph}("*")
    qtree = QueryAtom()
    anyatom = atomsymbol!(state, qtree)
    @test qtree == QueryAtom(Tuple{Int,Int}[], [qanytrue()]) # anything matches
end

@testset "atomprop" begin
    # inside of []
    state = SMARTSParser{SMARTSMolGraph}("s")
    qtree = QueryAtom()
    aroms = atomprop!(state, qtree)
    @test aroms == 1
    @test qtree == QueryAtom(
        [(1, 2), (1, 3)], [qand(), qeq(:symbol, "S"), qtrue(:isaromatic)])

    state = SMARTSParser{SMARTSMolGraph}("123")
    qtree = QueryAtom()
    iso = atomprop!(state, qtree)
    @test iso == 1
    @test qtree == QueryAtom(Tuple{Int,Int}[], [qeq(:mass, "123")])

    state = SMARTSParser{SMARTSMolGraph}("H41")
    qtree = QueryAtom()
    h4 = atomprop!(state, qtree)
    @test h4 == 1
    @test qtree == QueryAtom(Tuple{Int,Int}[], [qeq(:total_hydrogens, "4")])
    @test state.pos == 3  # Next, read '1'

    state = SMARTSParser{SMARTSMolGraph}("X2")
    qtree = QueryAtom()
    x2 = atomprop!(state, qtree)
    @test x2 == 1
    @test qtree == QueryAtom(Tuple{Int,Int}[], [qeq(:connectivity, "2")])

    state = SMARTSParser{SMARTSMolGraph}("Xe")  # not connectivity
    qtree = QueryAtom()
    xe = atomprop!(state, qtree)
    @test xe == 1
    @test qtree == QueryAtom(Tuple{Int,Int}[], [qeq(:symbol, "Xe")])

    state = SMARTSParser{SMARTSMolGraph}("Na")  # not nitrogen
    qtree = QueryAtom()
    na = atomprop!(state, qtree)
    @test na == 1
    @test qtree == QueryAtom(Tuple{Int,Int}[], [qeq(:symbol, "Na")])

    state = SMARTSParser{SMARTSMolGraph}("Yv2")
    qtree = QueryAtom()
    yatom = atomprop!(state, qtree)
    @test yatom == 1
    @test qtree == QueryAtom(Tuple{Int,Int}[], [qeq(:symbol, "Y")])
    yval = atomprop!(state, qtree)
    @test yval == 2
    qtree.vprops[2] == qeq(:valence, "2")

    state = SMARTSParser{SMARTSMolGraph}("+23")
    qtree = QueryAtom()
    chg1 = atomprop!(state, qtree)
    @test chg1 == 1
    @test qtree == QueryAtom(Tuple{Int,Int}[], [qeq(:charge, "2")])
    @test state.pos == 3  # Next, read '3'

    state = SMARTSParser{SMARTSMolGraph}("----+")
    qtree = QueryAtom()
    chg2 = atomprop!(state, qtree)
    @test chg2 == 1
    @test qtree == QueryAtom(Tuple{Int,Int}[], [qeq(:charge, "-4")])
    @test state.pos == 5  # Next, read '+'

    state = SMARTSParser{SMARTSMolGraph}("@+")
    qtree = QueryAtom()
    stereo1 = atomprop!(state, qtree)
    @test stereo1 == 1
    @test qtree == QueryAtom(Tuple{Int,Int}[], [qeq(:stereo, "anticlockwise")])
    @test state.pos == 2  # Next, read '+'

    state = SMARTSParser{SMARTSMolGraph}("@@?")  # clockwise or unspecified (racemate)
    qtree = QueryAtom()
    stereo2 = atomprop!(state, qtree)
    @test stereo2 == 1
    @test qtree == QueryAtom([(1, 2)], [qnot(), qeq(:stereo, "anticlockwise")])
    @test state.pos == 4

    state = SMARTSParser{SMARTSMolGraph}("\$([CH2]=*)")
    qtree = QueryAtom()
    rec = atomprop!(state, qtree)
    @test rec == 1
    @test qtree == QueryAtom(Tuple{Int,Int}[], [qeq(:recursive, "[CH2]=*")])
    @test state.pos == 11

    state = SMARTSParser{SMARTSMolGraph}("Uup")  # not supported yet
    qtree = QueryAtom()
    uup = atomprop!(state, qtree)
    @test uup == 1
    @test_broken qtree == QueryAtom(Tuple{Int,Int}[], [qeq(:symbol, "Uup")])
end

@testset "smarts" begin
    state = SMARTSParser{SMARTSMolGraph}("")
    null = atom!(state)
    @test isempty(null)

    state = SMARTSParser{SMARTSMolGraph}("a")
    anyarom = atom!(state)
    @test only(anyarom) isa QueryAtom
    @test only(anyarom) == QueryAtom(Tuple{Int,Int}[], [qtrue(:isaromatic)])

    state = SMARTSParser{SMARTSMolGraph}("[]")
    @test_throws ErrorException atom!(state)

    state = SMARTSParser{SMARTSMolGraph}("[*]")
    anyatom = atom!(state)
    @test only(anyatom) == QueryAtom(Tuple{Int,Int}[], [qanytrue()])

    state = SMARTSParser{SMARTSMolGraph}("[#16]")
    no16 = atom!(state)
    @test only(no16) == QueryAtom(Tuple{Int,Int}[], [qeq(:symbol, "S")])

    state = SMARTSParser{SMARTSMolGraph}("[CH2]")
    methylene = atom!(state)
    @test only(methylene) == QueryAtom(
        [(1, 2), (1, 3), (2, 4), (2, 5), (5, 6)],
        [qand(), qand(), qeq(:total_hydrogens, "2"), qeq(:symbol, "C"), qnot(), qtrue(:isaromatic)])

    state = SMARTSParser{SMARTSMolGraph}("[!C;R]")
    ringnotalp = atom!(state)
    @test only(ringnotalp) == QueryAtom(
        [(1, 2), (1, 3), (2, 4), (3, 5), (5, 6), (5, 7), (7, 8)],
        [qand(), qnot(), qnot(), qeq(:ring_count, "0"), qand(),
        qeq(:symbol, "C"), qnot(), qtrue(:isaromatic)])

    state = SMARTSParser{SMARTSMolGraph}("[n&H1]")
    nh1 = atom!(state)
    @test only(nh1) == QueryAtom(
        [(1, 2), (1, 3), (2, 4), (2, 5)],
        [qand(), qand(), qeq(:total_hydrogens, "1"),
        qeq(:symbol, "N"), qtrue(:isaromatic)])

    state = SMARTSParser{SMARTSMolGraph}("[*r6]")
    sixmem = atom!(state)
    @test only(sixmem) == QueryAtom(
        [(1, 2), (1, 3)], [qand(), qanytrue(), qeq(:smallest_ring, "6")])

    state = SMARTSParser{SMARTSMolGraph}("[35*]")
    any35 = atom!(state)
    @test only(any35) == QueryAtom(
        [(1, 2), (1, 3)], [qand(), qanytrue(), qeq(:mass, "35")])

    state = SMARTSParser{SMARTSMolGraph}("[F,Cl,Br,I]")
    fourhalo = atom!(state)
    @test only(fourhalo) == QueryAtom(
        [(1, 2), (1, 3), (1, 4), (1, 5)],
        [qor(), qeq(:symbol, "F") ,qeq(:symbol, "Cl") ,qeq(:symbol, "Br") ,qeq(:symbol, "I")])
end

@testset "smiles" begin
    # outside of []
    state = SMILESParser{SMILESMolGraph}("Br")
    br = atom!(state)
    @test only(br) == SMILESAtom(;symbol=:Br)

    state = SMILESParser{SMILESMolGraph}("se")  # only Se inside [] can be recognized
    se = atom!(state)
    @test only(se) == SMILESAtom(;symbol=:S, isaromatic=true)
    @test state.pos == 2  # Next, read 'e'

    state = SMILESParser{SMILESMolGraph}("Na")  # Invalid SMILES (SMARTS take it as N-a)
    invalid_na = atom!(state)
    @test only(invalid_na) == SMILESAtom(;symbol=:N)

    # inside of []
    state = SMILESParser{SMILESMolGraph}("[Na]")
    valid_na = atom!(state)
    @test only(valid_na) == SMILESAtom(;symbol=:Na)

    state = SMILESParser{SMILESMolGraph}("[2H]")
    deu = atom!(state)
    @test only(deu) == SMILESAtom(;symbol=:H, mass=2)

    state = SMILESParser{SMILESMolGraph}("[H2]")
    hmol = atom!(state)
    @test hmol[1] == hmol[2]

    state = SMILESParser{SMILESMolGraph}("[H+]")
    proton = atom!(state)
    @test only(proton) == SMILESAtom(;symbol=:H, charge=1)

    state = SMILESParser{SMILESMolGraph}("[14c@@H]")
    iso = atom!(state)
    @test iso[1] == SMILESAtom(;mass=14, isaromatic=true, stereo=:clockwise)
    @test iso[2] == SMILESAtom(;symbol=:H)
    @test state.pos == 9

    state = SMILESParser{SMILESMolGraph}("[NH4+]")
    ammonium = atom!(state)
    @test ammonium[1] == SMILESAtom(;symbol=:N, charge=1)
    @test ammonium[5] == SMILESAtom(;symbol=:H)

    state = SMILESParser{SMILESMolGraph}("[bH+]")
    aromb = atom!(state)
    @test aromb[1] == SMILESAtom(;symbol=:B, isaromatic=true, charge=1)
    @test aromb[2] == SMILESAtom(;symbol=:H)

    state = SMILESParser{SMILESMolGraph}("[as]")
    aromas = atom!(state)
    @test only(aromas) == SMILESAtom(;symbol=:As, isaromatic=true)

    state = SMILESParser{SMILESMolGraph}("[Zn++]")
    zn = atom!(state)
    @test only(zn) == SMILESAtom(;symbol=:Zn, charge=2)

    state = SMILESParser{SMILESMolGraph}("[O-2]")
    ox = atom!(state)
    @test only(ox) == SMILESAtom(;symbol=:O, charge=-2)
end

end # smiles.atom
