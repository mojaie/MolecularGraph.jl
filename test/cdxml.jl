
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
    @test atom_num_hs(mol) == [0]
    @test valence(mol) == [0]

    sdf = printv2mol(mol)
    @test contains(sdf, "C   0  0  0  0  0 15")  # vvv 15: zero valence

    mol2 = sdftomol(IOBuffer(sdf))  # vvv is read by SDFAtom (ignored by CDXMLAtom)
    @test atom_valence(mol2) == [0]
    @test round(Int, standard_weight(mol2)) == 12
end

@testset "NumHydrogens and valence" begin
    # methylene (carbene) label "CH2" bonded to two methyls, and a methyl with NumHydrogens="3"
    cdxml_string = """
    <?xml version="1.0" encoding="UTF-8" ?>
    <CDXML><page><fragment id="1">
    <n id="2" p="0 0"/>
    <n id="3" p="10 0" NumHydrogens="0"/>
    <n id="4" p="20 0" NumHydrogens="3"/>
    <b id="5" B="2" E="3"/>
    <b id="6" B="3" E="4"/>
    </fragment></page></CDXML>
    """
    mol = cdxmltomol(IOBuffer(cdxml_string))
    @test atom_num_hs(mol) == [-1, 0, 3]
    @test valence(mol) == [4, 2, 4]
    @test lone_pair(mol) == [0, 1, 0]
    @test molecular_formula(mol) == "C3H6"

    # non-organic elements
    sn_string = """
    <?xml version="1.0" encoding="UTF-8" ?>
    <CDXML><page><fragment id="1">
    <n id="2" p="0 0" Element="50" NumHydrogens="0"/>
    <n id="3" p="10 0" Element="50" NumHydrogens="2"/>
    <n id="4" p="20 0"/>
    <b id="5" B="3" E="4"/>
    </fragment></page></CDXML>
    """
    snmol = cdxmltomol(IOBuffer(sn_string))
    @test atom_num_hs(snmol) == [0, 2, -1]
    @test valence(snmol) == [0, 3, 4]
    @test snmol[:logs]["warning_valence"] == ""

    # aromatic atoms with NumHydrogens (pyrrole)
    pyrrole_string = """
    <?xml version="1.0" encoding="UTF-8" ?>
    <CDXML><page><fragment id="1">
    <n id="2" p="0 0" Element="7" NumHydrogens="1"/>
    <n id="3" p="10 0"/>
    <n id="4" p="20 0"/>
    <n id="5" p="20 10"/>
    <n id="6" p="10 10"/>
    <b id="7" B="2" E="3" Order="1.5"/>
    <b id="8" B="3" E="4" Order="1.5"/>
    <b id="9" B="4" E="5" Order="1.5"/>
    <b id="10" B="5" E="6" Order="1.5"/>
    <b id="11" B="6" E="2" Order="1.5"/>
    </fragment></page></CDXML>
    """
    pyrrole = cdxmltomol(IOBuffer(pyrrole_string))
    @test sum(bond_order(pyrrole) .== 2) == 2
    @test valence(pyrrole) == [3, 4, 4, 4, 4]
    @test molecular_formula(pyrrole) == "C4H5N"
end

@testset "hydrogen counts in molfiles" begin
    config = Dict{Symbol,Any}(
        :on_init => MolecularGraph.cdxml_on_init!, :on_update => MolecularGraph.cdxml_on_update!)
    readv3(atomline) = sdftomol(CDXMLMolGraph, IOBuffer("""

      test

      0  0  0     0  0            999 V3000
    M  V30 BEGIN CTAB
    M  V30 COUNTS 2 1 0 0 0
    M  V30 BEGIN ATOM
    M  V30 1 C 0.0 0.0 0.0 0
    $(atomline)
    M  V30 END ATOM
    M  V30 BEGIN BOND
    M  V30 1 1 1 2
    M  V30 END BOND
    M  V30 END CTAB
    M  END
    """); config)
    # HCOUNT: 0 = automatic, -1 = no hydrogens, n > 0 = n hydrogens
    @test atom_num_hs(readv3("M  V30 2 N 1.0 0.0 0.0 0")) == [-1, -1]
    @test atom_num_hs(readv3("M  V30 2 N 1.0 0.0 0.0 0 HCOUNT=-1")) == [-1, 0]
    @test valence(readv3("M  V30 2 N 1.0 0.0 0.0 0 HCOUNT=-1")) == [4, 1]
    @test atom_num_hs(readv3("M  V30 2 N 1.0 0.0 0.0 0 HCOUNT=2")) == [-1, 2]
    nh3 = readv3("M  V30 2 N 1.0 0.0 0.0 0 CHG=1 HCOUNT=3")
    @test atom_num_hs(nh3) == [-1, 3]
    @test implicit_hydrogens(nh3) == [3, 3]

    # V2000 hhh: n + 1 hydrogens, ignored by the default SDFAtom
    v2 = """

      test

      2  1  0  0  0  0  0  0  0  0999 V2000
        0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        1.0000    0.0000    0.0000 N   0  0  0  2  0  0  0  0  0  0  0  0
      1  2  1  0  0  0  0
    M  END
    """
    nh = sdftomol(CDXMLMolGraph, IOBuffer(v2); config)
    @test atom_num_hs(nh) == [-1, 1]
    @test implicit_hydrogens(nh) == [3, 1]
    @test implicit_hydrogens(sdftomol(IOBuffer(v2))) == [3, 2]
    # the property block (M  CHG) keeps the number of hydrogens from the atom block
    v2chg = replace(v2, "M  END" => "M  CHG  1   2   1\nM  END")
    nhchg = sdftomol(CDXMLMolGraph, IOBuffer(v2chg); config)
    @test atom_charge(nhchg) == [0, 1]
    @test atom_num_hs(nhchg) == [-1, 1]
    @test implicit_hydrogens(nhchg) == [3, 1]
end

end # cdxmlreader
