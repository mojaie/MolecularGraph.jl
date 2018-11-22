
@testset "smarts.atom" begin

@testset "atomsymbol" begin
    state = SmartsParserState("Cl")
    Chloride = atomsymbol!(state)
    @test Chloride == (:Cl, false)
    @test state.pos == 3

    state = SmartsParserState("Cr")
    Chromium = atomsymbol!(state)
    @test Chromium == (:C, false)
    @test state.pos == 2

    state = SmartsParserState("p")
    aromp = atomsymbol!(state)
    @test aromp == (:P, true)
end

@testset "atomprop" begin
    state = SmartsParserState("s")
    aroms = atomprop!(state)
    @test aroms == (:and => (:Symbol => :S, :Aromatic => true))

    state = SmartsParserState("123")
    iso = atomprop!(state)
    @test iso == (:Mass => 123)

    state = SmartsParserState("2H]")
    forward!(state)
    hyd = atomprop!(state)
    @test hyd == (:Symbol => :H)

    state = SmartsParserState("!H]")
    forward!(state)
    noth = atomprop!(state)
    @test noth == (:H_Count => 1)

    state = SmartsParserState("H+")
    proton = atomprop!(state)
    @test proton == (:and => (:Symbol => :H, :Charge => 1))

    state = SmartsParserState("H41")
    H4 = atomprop!(state)
    @test H4 == (:H_Count => 4)
    @test state.pos == 3

    state = SmartsParserState("X2")
    X2 = atomprop!(state)
    @test X2 == (:Connectivity => 2)

    state = SmartsParserState("Xe")
    Xe = atomprop!(state)
    @test Xe == (:Symbol => :Xe)
    @test state.pos == 3

    state = SmartsParserState("Na")
    Na = atomprop!(state)
    @test Na == (:Symbol => :Na)
    @test state.pos == 3

    state = SmartsParserState("Yv2")
    Yval = atomprop!(state)
    @test Yval == (:Symbol => :Y)

    state = SmartsParserState("+23")
    chg1 = atomprop!(state)
    @test chg1 == (:Charge => 2)
    @test state.pos == 3

    state = SmartsParserState("----+")
    chg2 = atomprop!(state)
    @test chg2 == (:Charge => -4)
    @test state.pos == 5

    state = SmartsParserState("@+")
    stereo1 = atomprop!(state)
    @test stereo1 == (:stereo => 1)
    @test state.pos == 2

    state = SmartsParserState("@@?")
    stereo4 = atomprop!(state)
    @test stereo4 == (:stereo => 4)
    @test state.pos == 4

    state = SmartsParserState("\$([CH2]=*)")
    rec = atomprop!(state)
    @test rec == (:recursive => "[CH2]=*")
    @test state.pos == 11
end

@testset "atom" begin
    state = SmilesParserState("Br")
    br = atom!(state)
    @test br.symbol == :Br

    state = SmilesParserState("[14c@@H]")
    iso = atom!(state)
    @test iso.symbol == :C
    @test iso.mass == 14
    @test iso.isaromatic
    @test iso.stereo == 2
    @test state.pos == 9

    state = SmilesParserState("[Zn++]")
    zn = atom!(state)
    @test zn.symbol == :Zn
    @test zn.charge == 2

    state = SmilesParserState("[O-2]")
    ox = atom!(state)
    @test ox.symbol == :O
    @test ox.charge == -2
end

@testset "smartsatom" begin
    state = SmartsParserState("c")
    aromc = atom!(state)
    @test aromc.query == (:and => (:Symbol => :C, :Aromatic => true))

    state = SmartsParserState("[#16]")
    no16 = atom!(state)
    @test no16.query == (:Symbol => :S)
end

end # smiles.atom
