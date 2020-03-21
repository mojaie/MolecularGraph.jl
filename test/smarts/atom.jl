
using MolecularGraph: atomsymbol!, atomprop!, atom!

@testset "smarts.atom" begin

@testset "atomsymbol" begin
    state = SmartsParser("Cl", false)
    Chloride = atomsymbol!(state)
    @test Chloride == (:and => (:atomsymbol => :Cl, :isaromatic => false))
    @test state.pos == 3

    state = SmartsParser("Cr", false)
    Chromium = atomsymbol!(state)
    @test Chromium == (:and => (:atomsymbol => :C, :isaromatic => false))
    @test state.pos == 2

    state = SmartsParser("p", false)
    aromp = atomsymbol!(state)
    @test aromp == (:and => (:atomsymbol => :P, :isaromatic => true))
end

@testset "atomprop" begin
    state = SmartsParser("s", false)
    aroms = atomprop!(state)
    @test aroms == (:and => (:atomsymbol => :S, :isaromatic => true))

    state = SmartsParser("123", false)
    iso = atomprop!(state)
    @test iso == (:mass => 123)

    state = SmartsParser("2H]", false)
    forward!(state)
    hyd = atomprop!(state)
    @test hyd == (:atomsymbol => :H)

    state = SmartsParser("!H]", false)
    forward!(state)
    noth = atomprop!(state)
    @test noth == (:hcount => 1)

    state = SmartsParser("H+", false)
    proton = atomprop!(state)
    @test proton == (:and => (:atomsymbol => :H, :charge => 1))

    state = SmartsParser("H41", false)
    H4 = atomprop!(state)
    @test H4 == (:hcount => 4)
    @test state.pos == 3

    state = SmartsParser("X2", false)
    X2 = atomprop!(state)
    @test X2 == (:connectivity => 2)

    state = SmartsParser("Xe", false)
    Xe = atomprop!(state)
    @test Xe == (:atomsymbol => :Xe)
    @test state.pos == 3

    state = SmartsParser("Na", false)
    Na = atomprop!(state)
    @test Na == (:atomsymbol => :Na)
    @test state.pos == 3

    state = SmartsParser("Yv2", false)
    Yval = atomprop!(state)
    @test Yval == (:atomsymbol => :Y)

    state = SmartsParser("+23", false)
    chg1 = atomprop!(state)
    @test chg1 == (:charge => 2)
    @test state.pos == 3

    state = SmartsParser("----+", false)
    chg2 = atomprop!(state)
    @test chg2 == (:charge => -4)
    @test state.pos == 5

    state = SmartsParser("@+", false)
    stereo1 = atomprop!(state)
    @test stereo1 == (:stereo => :anticlockwise)
    @test state.pos == 2

    state = SmartsParser("@@?", false)
    stereo4 = atomprop!(state)
    @test stereo4 == (:not => (:stereo => :anticlockwise))
    @test state.pos == 4

    state = SmartsParser("\$([CH2]=*)", false)
    rec = atomprop!(state)
    @test rec == (:recursive => "[CH2]=*")
    @test state.pos == 11
end

@testset "atom" begin
    state = SmilesParser("Br", true)
    br = atom!(state)
    @test br.symbol == :Br

    state = SmilesParser("[14c@@H]", true)
    iso = atom!(state)
    @test iso.symbol == :C
    @test iso.mass == 14
    @test iso.isaromatic
    @test iso.stereo == :clockwise
    @test state.pos == 9

    state = SmilesParser("[Zn++]", true)
    zn = atom!(state)
    @test zn.symbol == :Zn
    @test zn.charge == 2

    state = SmilesParser("[O-2]", true)
    ox = atom!(state)
    @test ox.symbol == :O
    @test ox.charge == -2
end

@testset "smartsatom" begin
    state = SmartsParser("c", false)
    aromc = atom!(state)
    @test aromc.query == (:and => (:atomsymbol => :C, :isaromatic => true))

    state = SmartsParser("[#16]", false)
    no16 = atom!(state)
    @test no16.query == (:atomsymbol => :S)
end

end # smiles.atom
