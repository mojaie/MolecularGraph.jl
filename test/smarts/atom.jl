
@testset "smarts.atom" begin

@testset "atomsymbol" begin
    state = SmartsParser{SMARTS}("Cl", false)
    Chloride = atomsymbol!(state)
    @test Chloride == (:and => (:Symbol => :Cl, :Aromatic => false))
    @test state.pos == 3

    state = SmartsParser{SMARTS}("Cr", false)
    Chromium = atomsymbol!(state)
    @test Chromium == (:and => (:Symbol => :C, :Aromatic => false))
    @test state.pos == 2

    state = SmartsParser{SMARTS}("p", false)
    aromp = atomsymbol!(state)
    @test aromp == (:and => (:Symbol => :P, :Aromatic => true))
end

@testset "atomprop" begin
    state = SmartsParser{SMARTS}("s", false)
    aroms = atomprop!(state)
    @test aroms == (:and => (:Symbol => :S, :Aromatic => true))

    state = SmartsParser{SMARTS}("123", false)
    iso = atomprop!(state)
    @test iso == (:Mass => 123)

    state = SmartsParser{SMARTS}("2H]", false)
    forward!(state)
    hyd = atomprop!(state)
    @test hyd == (:Symbol => :H)

    state = SmartsParser{SMARTS}("!H]", false)
    forward!(state)
    noth = atomprop!(state)
    @test noth == (:H_Count => 1)

    state = SmartsParser{SMARTS}("H+", false)
    proton = atomprop!(state)
    @test proton == (:and => (:Symbol => :H, :Charge => 1))

    state = SmartsParser{SMARTS}("H41", false)
    H4 = atomprop!(state)
    @test H4 == (:H_Count => 4)
    @test state.pos == 3

    state = SmartsParser{SMARTS}("X2", false)
    X2 = atomprop!(state)
    @test X2 == (:Connectivity => 2)

    state = SmartsParser{SMARTS}("Xe", false)
    Xe = atomprop!(state)
    @test Xe == (:Symbol => :Xe)
    @test state.pos == 3

    state = SmartsParser{SMARTS}("Na", false)
    Na = atomprop!(state)
    @test Na == (:Symbol => :Na)
    @test state.pos == 3

    state = SmartsParser{SMARTS}("Yv2", false)
    Yval = atomprop!(state)
    @test Yval == (:Symbol => :Y)

    state = SmartsParser{SMARTS}("+23", false)
    chg1 = atomprop!(state)
    @test chg1 == (:Charge => 2)
    @test state.pos == 3

    state = SmartsParser{SMARTS}("----+", false)
    chg2 = atomprop!(state)
    @test chg2 == (:Charge => -4)
    @test state.pos == 5

    state = SmartsParser{SMARTS}("@+", false)
    stereo1 = atomprop!(state)
    @test stereo1 == (:stereo => 1)
    @test state.pos == 2

    state = SmartsParser{SMARTS}("@@?", false)
    stereo4 = atomprop!(state)
    @test stereo4 == (:stereo => 4)
    @test state.pos == 4

    state = SmartsParser{SMARTS}("\$([CH2]=*)", false)
    rec = atomprop!(state)
    @test rec == (:recursive => "[CH2]=*")
    @test state.pos == 11
end

@testset "atom" begin
    state = SmartsParser{SMILES}("Br", true)
    br = atom!(state)
    @test br.symbol == :Br

    state = SmartsParser{SMILES}("[14c@@H]", true)
    iso = atom!(state)
    @test iso.symbol == :C
    @test iso.mass == 14
    @test iso.isaromatic
    @test iso.stereo == 2
    @test state.pos == 9

    state = SmartsParser{SMILES}("[Zn++]", true)
    zn = atom!(state)
    @test zn.symbol == :Zn
    @test zn.charge == 2

    state = SmartsParser{SMILES}("[O-2]", true)
    ox = atom!(state)
    @test ox.symbol == :O
    @test ox.charge == -2
end

@testset "smartsatom" begin
    state = SmartsParser{SMARTS}("c", false)
    aromc = atom!(state)
    @test aromc.query == (:and => (:Symbol => :C, :Aromatic => true))

    state = SmartsParser{SMARTS}("[#16]", false)
    no16 = atom!(state)
    @test no16.query == (:Symbol => :S)
end

end # smiles.atom
