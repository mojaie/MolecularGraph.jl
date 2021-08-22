
using MolecularGraph: atomsymbol!, atomprop!, atom!

@testset "smarts.atom" begin

@testset "atomsymbol" begin
    # outside of []
    state = SmartsParser("Cl", false)
    Chloride = atomsymbol!(state)
    @test Chloride == QueryFormula(:and, Set([
        QueryFormula(:atomsymbol, :Cl),
        QueryFormula(:isaromatic, false)
    ]))
    @test state.pos == 3

    state = SmartsParser("Cr", false)
    Chromium = atomsymbol!(state)
    @test Chromium == QueryFormula(:and, Set([
        QueryFormula(:atomsymbol, :C),
        QueryFormula(:isaromatic, false)
    ]))
    @test state.pos == 2

    state = SmartsParser("p", false)
    aromp = atomsymbol!(state)
    
    @test aromp == QueryFormula(:and, Set([
        QueryFormula(:atomsymbol, :P),
        QueryFormula(:isaromatic, true)
    ]))

    state = SmartsParser("*", false)
    anyatom = atomsymbol!(state)
    @test anyatom == QueryFormula(:any, true)
end

@testset "atomprop" begin
    # inside of []
    state = SmartsParser("s", false)
    aroms = atomprop!(state)
    @test aroms == QueryFormula(:and, Set([
        QueryFormula(:atomsymbol, :S),
        QueryFormula(:isaromatic, true)
    ]))

    state = SmartsParser("123", false)
    iso = atomprop!(state)
    @test iso == QueryFormula(:mass, 123)

    state = SmartsParser("H41", false)
    H4 = atomprop!(state)
    @test H4 == QueryFormula(:hydrogenconnected, 4)
    @test state.pos == 3

    state = SmartsParser("X2", false)
    X2 = atomprop!(state)
    @test X2 == QueryFormula(:connectivity, 2)

    state = SmartsParser("Xe", false)
    Xe = atomprop!(state)
    @test Xe == QueryFormula(:atomsymbol, :Xe)
    @test state.pos == 3

    state = SmartsParser("Na", false)
    Na = atomprop!(state)
    @test Na == QueryFormula(:atomsymbol, :Na)
    @test state.pos == 3

    state = SmartsParser("Yv2", false)
    Yval = atomprop!(state)
    @test Yval == QueryFormula(:atomsymbol, :Y)

    state = SmartsParser("+23", false)
    chg1 = atomprop!(state)
    @test chg1 == QueryFormula(:charge, 2)
    @test state.pos == 3

    state = SmartsParser("----+", false)
    chg2 = atomprop!(state)
    @test chg2 == QueryFormula(:charge, -4)
    @test state.pos == 5

    state = SmartsParser("@+", false)
    stereo1 = atomprop!(state)
    @test stereo1 == QueryFormula(:stereo, :anticlockwise)
    @test state.pos == 2

    state = SmartsParser("@@?", false)
    stereo4 = atomprop!(state)
    @test stereo4 == QueryFormula(:not, QueryFormula(:stereo, :anticlockwise))
    @test state.pos == 4

    state = SmartsParser("\$([CH2]=*)", false)
    rec = atomprop!(state)
    @test rec == QueryFormula(:recursive, "[CH2]=*")
    @test state.pos == 11
end

@testset "atom" begin
    state = SmilesParser("Br", false)
    br = atom!(state)[1]
    @test br.symbol == :Br

    state = SmilesParser("[2H]", false)
    hyd = atom!(state)
    @test length(hyd) == 1
    @test hyd[1].mass == 2

    state = SmilesParser("[H2]", false)
    hmol = atom!(state)
    @test length(hmol) == 2

    state = SmilesParser("[H+]", false)
    proton = atom!(state)[1]
    @test proton.charge == 1

    state = SmilesParser("[14c@@H]", false)
    iso = atom!(state)
    @test iso[1].symbol == :C
    @test iso[1].mass == 14
    @test iso[1].isaromatic
    @test iso[1].stereo == :clockwise
    @test iso[2].symbol == :H
    @test length(iso) == 2
    @test state.pos == 9

    state = SmilesParser("[NH4+]", false)
    ammonium = atom!(state)
    @test length(ammonium) == 5
    @test ammonium[1].charge == 1
    @test ammonium[5].symbol == :H

    state = SmilesParser("[bH+]", false)
    aromB = atom!(state)
    @test aromB[1].symbol == :B
    @test aromB[1].isaromatic
    @test aromB[1].charge == 1
    @test aromB[2].symbol == :H

    state = SmilesParser("[as]", false)
    aromAs = atom!(state)[1]
    @test aromAs.symbol == :As
    @test aromAs.isaromatic

    state = SmilesParser("[Zn++]", false)
    zn = atom!(state)[1]
    @test zn.symbol == :Zn
    @test zn.charge == 2

    state = SmilesParser("[O-2]", false)
    ox = atom!(state)[1]
    @test ox.symbol == :O
    @test ox.charge == -2

end

@testset "smartsatom" begin
    state = SmartsParser("c", false)
    aromc = atom!(state)[1]
    @test aromc.query == QueryFormula(:and, Set([
        QueryFormula(:atomsymbol, :C),
        QueryFormula(:isaromatic, true)
    ]))

    state = SmartsParser("a", false)
    anyarom = atom!(state)[1]
    @test anyarom.query == QueryFormula(:isaromatic, true)

    state = SmartsParser("[]", false)
    @test_throws AssertionError atom!(state)

    state = SmartsParser("[*]", false)
    anyatom = atom!(state)[1]
    @test anyatom.query == QueryFormula(:any, true)

    state = SmartsParser("[#16]", false)
    no16 = atom!(state)[1]
    @test no16.query == QueryFormula(:atomsymbol, :S)

    state = SmartsParser("[CH2]", false)
    methylene = atom!(state)[1]
    @test methylene.query == QueryFormula(:and, Set([
        QueryFormula(:atomsymbol, :C),
        QueryFormula(:isaromatic, false),
        QueryFormula(:hydrogenconnected, 2)
    ]))

    state = SmartsParser("[!C;R]", false)
    ringnotalp = atom!(state)[1]
    @test ringnotalp.query == QueryFormula(:and, Set([
        QueryFormula(:or, Set([
            QueryFormula(:not, QueryFormula(:atomsymbol, :C)),
            QueryFormula(:isaromatic, true)
        ])),
        QueryFormula(:not, QueryFormula(:sssrcount, 0))
    ]))

    state = SmartsParser("[n&H1]", false)
    nh1 = atom!(state)[1]
    @test nh1.query == QueryFormula(:and, Set([
        QueryFormula(:atomsymbol, :N),
        QueryFormula(:isaromatic, true),
        QueryFormula(:hydrogenconnected, 1)
    ]))

    state = SmartsParser("[*r6]", false)
    sixmem = atom!(state)[1]
    @test sixmem.query == QueryFormula(:sssrsizes, 6)

    state = SmartsParser("[35*]", false)
    any35 = atom!(state)[1]
    @test any35.query == QueryFormula(:mass, 35)

    state = SmartsParser("[F,Cl,Br,I]", false)
    fourhalo = atom!(state)[1]
    @test fourhalo.query == QueryFormula(:or, Set([
        QueryFormula(:atomsymbol, :F),
        QueryFormula(:atomsymbol, :Cl),
        QueryFormula(:atomsymbol, :Br),
        QueryFormula(:atomsymbol, :I)
    ]))
end

end # smiles.atom
