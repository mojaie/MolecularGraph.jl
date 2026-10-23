
@testset "cdxmlreader" begin

@testset "aromatic ring" begin
    # cdxml codes aromaticity via bond order 1.5
    cdxml_string = """
    <?xml version="1.0" encoding="UTF-8"?>
    <CDXML>
        <page>
        <fragment id="1">
            <n id="2" Element="6" p="100 100"/>
            <n id="3" Element="6" p="120 110"/>
            <n id="4" Element="6" p="120 130"/>
            <n id="5" Element="6" p="100 140"/>
            <n id="6" Element="6" p="80 130"/>
            <n id="7" Element="6" p="80 110"/>
            <b id="8" B="2" E="3" Order="1.5"/>
            <b id="9" B="3" E="4" Order="1.5"/>
            <b id="10" B="4" E="5" Order="1.5"/>
            <b id="11" B="5" E="6" Order="1.5"/>
            <b id="12" B="6" E="7" Order="1.5"/>
            <b id="13" B="7" E="2" Order="1.5"/>
        </fragment>
        </page>
    </CDXML>
    """

    mol = cdxmltomol(IOBuffer(cdxml_string))
    @test atom_symbol(mol) == fill(:C, 6)
    @test sum(bond_order(mol)) == 6 * 1.5
    @test mol[1].coords == [100, 100, 0]

    benzene_cdxml_aromatic = mol
    cdxml_string = """
    <?xml version="1.0" encoding="UTF-8"?>
    <CDXML>
        <page>
        <fragment id="1">
            <n id="2" Element="6" p="100 100"/>
            <n id="3" Element="6" p="120 110"/>
            <n id="4" Element="6" p="120 130"/>
            <n id="5" Element="6" p="100 140"/>
            <n id="6" Element="6" p="80 130"/>
            <n id="7" Element="6" p="80 110"/>
            <b id="8" B="2" E="3" Order="1"/>
            <b id="9" B="3" E="4" Order="2"/>
            <b id="10" B="4" E="5" Order="1"/>
            <b id="11" B="5" E="6" Order="2"/>
            <b id="12" B="6" E="7" Order="1"/>
            <b id="13" B="7" E="2" Order="2"/>
        </fragment>
        </page>
    </CDXML>
    """
    benzene_cdxml_classic = cdxmltomol(IOBuffer(cdxml_string))
    @test has_exact_match(benzene_cdxml_aromatic, benzene_cdxml_classic)
    
    benzene_smiles = smilestomol("c1ccccc1")
    @test has_exact_match(benzene_cdxml_aromatic, benzene_smiles)
end


@testset "stereo bond" begin
    cdxml_string = """
    <?xml version="1.0" encoding="UTF-8"?>
    <CDXML>
    <page>
        <fragment id="1">
        <n id="2" Element="6" p="100 100"/>
        <n id="3" Element="17" p="110 95"/>
        <n id="4" Element="35" p="110 105"/>
        <n id="5" Element="1" p="90 95"/>
        <n id="6" Element="1" p="90 105"/>
        <b id="7" B="2" E="3" Display="WedgeBegin"/>
        <b id="8" B="2" E="4" Display="Hash"/>
        <b id="9" B="2" E="5"/>
        <b id="10" B="2" E="6"/>
        </fragment>
    </page>
    </CDXML>
    """

    mol = cdxmltomol(IOBuffer(cdxml_string))
    @test [mol[e].notation for e in edges(mol)] == [1, 6, 0, 0]
    @test bond_order(mol) == fill(1, 4)
    @test atom_symbol(mol) == [:C, :Cl, :Br, :H, :H]
    @test mol[2].coords == [110, 95, 0]
end

@testset "NumHydrogens in charcoal" begin
    cdxml_string = """
    <?xml version="1.0" encoding="UTF-8" ?>
    <!DOCTYPE CDXML SYSTEM "https://static.chemistry.revvitycloud.com/cdxml/CDXML.dtd" >
    <CDXML><page><fragment
    id="3"
    BoundingBox="-3.60 -4.72 3.60 3.55"
    Z="1"
    ><n
    id="4"
    p="0 0"
    Z="2"
    Warning="An atom in this label has an invalid valence."
    NumHydrogens="0"
    NeedsClean="yes"
    AS="N"
    AtomNumber="0"
    AtomID="1"
    ><t
    p="-3.60 3.45"
    BoundingBox="-3.60 -4.72 3.60 3.55"
    LabelJustification="Left"
    ><s font="3" size="10" color="0" face="96">C</s></t></n></fragment></page></CDXML>
    """
    mol = cdxmltomol(IOBuffer(cdxml_string))

    @test mol isa CDXMLMolGraph
    @test round(Int, standard_weight(mol)) == 12
    
    sdf = printv2mol(mol)
    @test contains(sdf, "C   0  0  0  0  0 15")  # zero valence
    
    mol2 = sdftomol(MolecularGraph.CDXMLMolGraph, IOBuffer(sdf),
        config = Dict{Symbol,Any}(
            :on_init => MolecularGraph.cdxml_on_init!,
            :on_update => MolecularGraph.cdxml_on_update!
        )
    )

    @test mol isa CDXMLMolGraph
    @test round(Int, standard_weight(mol2)) == 12
end

@testset "fixed number of hydrogens in molfiles" begin
    cdxmlmol(molfile) = sdftomol(CDXMLMolGraph, IOBuffer(molfile);
        config=(on_init=MolecularGraph.cdxml_on_init!, on_update=MolecularGraph.cdxml_on_update!))
    molfile = "\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n" *
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0 15  0  0  0  0  0  0\nM  END\n"
    mol = cdxmlmol(molfile)
    @test mol[1].implicit_hs == 0
    @test implicit_hydrogens(mol) == [0]
    # V3000 with HCOUNT (no automatic hydrogens) and VAL
    ctab_atom_v3 = MolecularGraph.ctab_atom_v3
    @test ctab_atom_v3(CDXMLAtom, "M  V30 1 N 0 0 0 0 CHG=1 HCOUNT=2").implicit_hs == -1  # set after the bonds are read
    # round trip between the molecule types keeps the number of hydrogens
    atoms = "    0.0000    0.0000    0.0000 O   0  0  0  0  0  1  0  0  0  0  0  0\n"
    hydroxide = "\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n" * atoms * "M  CHG  1   1  -1\nM  END\n"
    ethyl = "\n  test\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n" *
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  3  0  0  0  0  0  0\n" *
        "    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\nM  END\n"
    for molfile in (molfile, hydroxide, ethyl)
        sdf = sdftomol(IOBuffer(molfile))
        cdxml = cdxmlmol(molfile)
        @test implicit_hydrogens(cdxml) == implicit_hydrogens(sdf)
        @test implicit_hydrogens(cdxmlmol(printv2mol(cdxml))) == implicit_hydrogens(sdf)
        @test implicit_hydrogens(sdftomol(IOBuffer(printv2mol(cdxml)))) == implicit_hydrogens(sdf)
        @test implicit_hydrogens(cdxmlmol(printv2mol(sdf))) == implicit_hydrogens(sdf)
    end
    for smiles in ("[C]", "[OH-]", "[CH2]C", "[O-]C")
        smol = smilestomol(smiles)
        remove_hydrogens!(smol)
        @test implicit_hydrogens(cdxmlmol(printv2mol(smol))) == implicit_hydrogens(smol)
    end
    @test implicit_hydrogens(sdftomol(IOBuffer(hydroxide))) == [1]
    @test implicit_hydrogens(sdftomol(IOBuffer(ethyl))) == [2, 3]
end

end # cdxmlreader
