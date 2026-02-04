
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
    @test [b.notation for (_, b) in mol.eprops] == [1, 6, 0, 0]
    @test bond_order(mol) == fill(1, 4)
    @test atom_symbol(mol) == [:C, :Cl, :Br, :H, :H]
    @test mol[2].coords == [110, 95, 0]
end

end # cdxmlreader
