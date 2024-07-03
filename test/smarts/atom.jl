
using MolecularGraph: atomsymbol!, atomprop!, atom!

@testset "smarts.atom" begin

@testset "atomsymbol" begin
    # outside of []
    state = SMARTSParser{SMARTSMolGraph}("Cl")
    cl = QueryTruthTable(atomsymbol!(state))
    @test cl == QueryTruthTable(v -> v[2] & ~v[1], [(:isaromatic,), (:symbol, :Cl)])
    @test state.pos == 3

    state = SMARTSParser{SMARTSMolGraph}("Cr")  # only Cr inside [] can be recognized
    cr = QueryTruthTable(atomsymbol!(state))
    @test cr == QueryTruthTable(v -> v[2] & ~v[1], [(:isaromatic,), (:symbol, :C)])
    @test state.pos == 2  # Next, read 'r'

    state = SMARTSParser{SMARTSMolGraph}("p")
    aromp = QueryTruthTable(atomsymbol!(state))
    @test aromp == QueryTruthTable(v -> v[2] & v[1], [(:isaromatic,), (:symbol, :P)])

    state = SMARTSParser{SMARTSMolGraph}("*")
    anyatom = QueryTruthTable(atomsymbol!(state))
    @test anyatom == QueryTruthTable(v -> true, [])  # anything matches
end

@testset "atomprop" begin
    # inside of []
    state = SMARTSParser{SMARTSMolGraph}("s")
    aroms = QueryTruthTable(atomprop!(state))
    @test aroms == QueryTruthTable(v -> v[2] & v[1], [(:isaromatic,), (:symbol, :S)])

    state = SMARTSParser{SMARTSMolGraph}("123")
    iso = QueryTruthTable(atomprop!(state))
    @test iso == QueryTruthTable(v -> v[1], [(:mass, 123)])

    state = SMARTSParser{SMARTSMolGraph}("H41")
    h4 = QueryTruthTable(atomprop!(state))
    @test h4 == QueryTruthTable(v -> v[1], [(:total_hydrogens, 4)])
    @test state.pos == 3  # Next, read '1'

    state = SMARTSParser{SMARTSMolGraph}("X2")
    x2 = QueryTruthTable(atomprop!(state))
    @test x2 == QueryTruthTable(v -> v[1], [(:connectivity, 2)])

    state = SMARTSParser{SMARTSMolGraph}("Xe")
    xe = QueryTruthTable(atomprop!(state))
    @test xe == QueryTruthTable(v -> v[1], [(:symbol, :Xe)])
    @test state.pos == 3

    state = SMARTSParser{SMARTSMolGraph}("Na")
    na = QueryTruthTable(atomprop!(state))
    @test na == QueryTruthTable(v -> v[1], [(:symbol, :Na)])
    @test state.pos == 3

    state = SMARTSParser{SMARTSMolGraph}("Yv2")
    yval = QueryTruthTable(atomprop!(state))
    @test yval == QueryTruthTable(v -> v[1], [(:symbol, :Y)])
    @test state.pos == 2  # Next, read 'v', '2' (valence: 2)

    state = SMARTSParser{SMARTSMolGraph}("+23")
    chg1 = QueryTruthTable(atomprop!(state))
    @test chg1 == QueryTruthTable(v -> v[1], [(:charge, 2)])
    @test state.pos == 3  # Next, read '3'

    state = SMARTSParser{SMARTSMolGraph}("----+")
    chg2 = QueryTruthTable(atomprop!(state))
    @test chg2 == QueryTruthTable(v -> v[1], [(:charge, -4)])
    @test state.pos == 5  # Next, read '+'

    state = SMARTSParser{SMARTSMolGraph}("@+")
    stereo1 = QueryTruthTable(atomprop!(state))
    @test stereo1 == QueryTruthTable(v -> v[1], [(:stereo, :anticlockwise)])
    @test state.pos == 2  # Next, read '+'

    state = SMARTSParser{SMARTSMolGraph}("@@?")  # clockwise or unspecified (racemate)
    stereo2 = QueryTruthTable(atomprop!(state))
    @test stereo2 == QueryTruthTable(v -> ~v[1], [(:stereo, :anticlockwise)])
    @test state.pos == 4

    state = SMARTSParser{SMARTSMolGraph}("\$([CH2]=*)")
    rec = QueryTruthTable(atomprop!(state))
    @test rec == QueryTruthTable(v -> v[1], [(:recursive, "[CH2]=*")])
    @test state.pos == 11
end

@testset "smilesatom" begin
    state = SMILESParser{SMILESMolGraph}("Br")
    br = atom!(state)[1]
    @test br[:symbol] === :Br

    state = SMILESParser{SMILESMolGraph}("Na")  # Invalid SMILES (SMARTS take it as N-a)
    invalid_na = atom!(state)[1]
    @test invalid_na[:symbol] === :N

    state = SMILESParser{SMILESMolGraph}("[Na]")
    valid_na = atom!(state)[1]
    @test valid_na[:symbol] === :Na

    state = SMILESParser{SMILESMolGraph}("[2H]")
    deu = atom!(state)
    @test length(deu) == 1
    @test deu[1][:mass] == 2

    state = SMILESParser{SMILESMolGraph}("[H2]")
    hmol = atom!(state)
    @test length(hmol) == 2

    state = SMILESParser{SMILESMolGraph}("[H+]")
    proton = atom!(state)[1]
    @test proton[:charge] == 1

    state = SMILESParser{SMILESMolGraph}("[14c@@H]")
    iso = atom!(state)
    @test iso[1][:symbol] === :C
    @test iso[1][:mass] == 14
    @test iso[1][:isaromatic]
    @test iso[1][:stereo] === :clockwise
    @test iso[2][:symbol] === :H
    @test length(iso) == 2
    @test state.pos == 9

    state = SMILESParser{SMILESMolGraph}("[NH4+]")
    ammonium = atom!(state)
    @test length(ammonium) == 5
    @test ammonium[1][:charge] == 1
    @test ammonium[5][:symbol] === :H

    state = SMILESParser{SMILESMolGraph}("[bH+]")
    aromb = atom!(state)
    @test aromb[1][:symbol] === :B
    @test aromb[1][:isaromatic]
    @test aromb[1][:charge] == 1
    @test aromb[2][:symbol] === :H

    state = SMILESParser{SMILESMolGraph}("[as]")
    aromas = atom!(state)[1]
    @test aromas[:symbol] === :As
    @test aromas[:isaromatic]

    state = SMILESParser{SMILESMolGraph}("[Zn++]")
    zn = atom!(state)[1]
    @test zn[:symbol] === :Zn
    @test zn[:charge] == 2

    state = SMILESParser{SMILESMolGraph}("[O-2]")
    ox = atom!(state)[1]
    @test ox[:symbol] === :O
    @test ox[:charge] == -2

end

@testset "smartsatom" begin
    SMARTSTT = MolGraph{Int,QueryTruthTable,QueryTruthTable}
    state = SMARTSParser{SMARTSTT}("c")
    aromc = atom!(state)[1]
    @test aromc == QueryTruthTable(v -> v[2] & v[1], [(:isaromatic,), (:symbol, :C)])

    state = SMARTSParser{SMARTSTT}("a")
    anyarom = atom!(state)[1]
    @test anyarom == QueryTruthTable(v -> v[1], [(:isaromatic,)])

    state = SMARTSParser{SMARTSTT}("[]")
    @test_throws ErrorException atom!(state)

    state = SMARTSParser{SMARTSTT}("[*]")
    anyatom = atom!(state)[1]
    @test anyatom == QueryTruthTable(v -> true, [])

    state = SMARTSParser{SMARTSTT}("[#16]")
    no16 = atom!(state)[1]
    @test no16 == QueryTruthTable(v -> v[1], [(:symbol, :S)])

    state = SMARTSParser{SMARTSTT}("[CH2]")
    methylene = atom!(state)[1]
    @test methylene == QueryTruthTable(
        v -> v[2] & ~v[1] & v[3], [(:isaromatic,), (:symbol, :C), (:total_hydrogens, 2)])

    state = SMARTSParser{SMARTSTT}("[!C;R]")
    ringnotalp = atom!(state)[1]
    @test ringnotalp == QueryTruthTable(
        v -> (~v[3] | v[1]) & ~v[2], [(:isaromatic,), (:ring_count, 0), (:symbol, :C)])

    state = SMARTSParser{SMARTSTT}("[n&H1]")
    nh1 = atom!(state)[1]
    @test nh1 == QueryTruthTable(
        v -> v[2] & v[1] & v[3], [(:isaromatic,), (:symbol, :N), (:total_hydrogens, 1)])

    state = SMARTSParser{SMARTSTT}("[*r6]")
    sixmem = atom!(state)[1]
    @test sixmem == QueryTruthTable(v -> v[1], [(:smallest_ring, 6)])

    state = SMARTSParser{SMARTSTT}("[35*]")
    any35 = atom!(state)[1]
    @test any35 == QueryTruthTable(v -> v[1], [(:mass, 35)])

    state = SMARTSParser{SMARTSTT}("[F,Cl,Br,I]")
    fourhalo = atom!(state)[1]
    @test fourhalo == QueryTruthTable(
        v -> v[1] | v[2] | v[3] | v[4],
        [(:symbol, :Br), (:symbol, :Cl), (:symbol, :F), (:symbol, :I)]
    )
end

end # smiles.atom
