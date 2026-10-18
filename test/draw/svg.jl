
@testset "svg" begin
    ASSETS_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets")
    TEST_DIR = joinpath(ASSETS_DIR, "test")
    # Expected SVG outputs. Set MOLGRAPH_UPDATE_SVG=true to regenerate them
    # after intentional changes in the drawing functions.
    EXPECTED_DIR = joinpath(TEST_DIR, "svg")
    UPDATE_SVG = get(ENV, "MOLGRAPH_UPDATE_SVG", "false") == "true"

    @testset "draw" begin
        demomol = sdftomol(open(joinpath(TEST_DIR, "demo.mol")))
        drawsvg(demomol)
        nullmol = sdftomol(open(joinpath(TEST_DIR, "null.mol")))
        @test count("<line", drawsvg(nullmol)) == 0
        # dest = open(joinpath(ASSETS_DIR, "image", "demo.svg"), "w")
        # write(dest, drawsvg(demomol))
    end

    @testset "structure" begin
        for file in ("aspirin_3d.sdf", "demo.mol", "biotin_2d.sdf", "DB00115.mol", "null.mol")
            mol = sdftomol(joinpath(TEST_DIR, file))
            svg = drawsvg(mol)
            # well-formed XML
            @test EzXML.nodename(EzXML.root(EzXML.parsexml(svg))) == "svg"
            # one text element for each visible atom label
            @test count("<text", svg) == count(is_atom_visible(mol))
        end
        biotin = sdftomol(joinpath(TEST_DIR, "biotin_2d.sdf"))
        @test count("<line", drawsvg(biotin)) == 53
    end

    @testset "options" begin
        mol = sdftomol(joinpath(TEST_DIR, "aspirin_3d.sdf"))
        default = drawsvg(mol)
        @test occursin("fill=\"rgb(255, 255, 255)\"", default)  # default bgcolor #FFF
        @test occursin("width=\"100%\" height=\"100%\" viewBox=", default)
        # background color
        svg = drawsvg(mol; bgcolor="#F0F0F0")
        @test occursin("fill=\"rgb(240, 240, 240)\"", svg)
        # atom highlight (default color #FDD835)
        svg = drawsvg(mol; atomhighlight=[1, 2])
        @test count("fill=\"rgb(253, 216, 53)\"", svg) == 2
        # bond highlight with a color name
        svg = drawsvg(mol; bondhighlight=[Edge(1, 5)], highlightcolor="red")
        @test count("stroke=\"rgb(255, 0, 0)\" stroke-width", svg) == 1
        @test count("stroke=\"rgb(255, 0, 0)\" stroke-width", drawsvg(mol; bondhighlight=[Edge(5, 1)], highlightcolor="red")) == 1  # reversed order
        @test_throws ArgumentError drawsvg(mol; bondhighlight=[Edge(1, 2)])  # no bond between 1 and 2
        # atom index
        svg = drawsvg(mol; atomindex=true)
        @test count("<g>", svg) == nv(mol)
        @test EzXML.nodename(EzXML.root(EzXML.parsexml(svg))) == "svg"
        # carbon labels
        nc = count(==(:C), atom_symbol(mol))
        @test count("<text", drawsvg(mol; show_carbon=:all)) == count("<text", default) + nc
        @test count("<g>", drawsvg(mol; show_carbon=:all, atomindex=true)) == nv(mol)
        # fixed size
        svg = drawsvg(mol; viewbox=false)
        @test !occursin("viewBox", svg)
        @test occursin(r"width=\"[0-9.]+\" height=\"[0-9.]+\"", svg)
    end

    @testset "parse_rgb" begin
        @test svgcolor(parse_rgb("#FFF")) == "rgb(255, 255, 255)"
        @test svgcolor(parse_rgb("#FDD835")) == "rgb(253, 216, 53)"
        @test svgcolor(parse_rgb("0xFF8800")) == "rgb(255, 136, 0)"
        @test svgcolor(parse_rgb(" SkyBlue ")) == "rgb(135, 206, 235)"
        @test svgcolor(parse_rgb(RGB(1, 0, 0))) == "rgb(255, 0, 0)"
        @test_throws ArgumentError parse_rgb("notacolor")
    end

    @testset "regression" begin
        # Molecules with 2D coordinates in the file (not depending on coordgen)
        cases = [
            ("demo.mol", "demo.svg", (;)),  # wedges, wavy bond, ring double bonds, charges
            ("biotin_2d.sdf", "biotin_2d.svg", (;)),  # stereocenters
            ("DB00115.mol", "DB00115.svg", (;)),  # large molecule with metal
            # 2D coordinates are the x and y of the 3D coordinates
            ("aspirin_3d.sdf", "aspirin_options.svg", (
                atomhighlight=[1, 2], bondhighlight=[Edge(1, 5)],
                atomindex=true, bgcolor="#F0F0F0"))
        ]
        UPDATE_SVG && mkpath(EXPECTED_DIR)
        for (molfile, svgfile, kwargs) in cases
            svg = drawsvg(sdftomol(joinpath(TEST_DIR, molfile)); kwargs...)
            path = joinpath(EXPECTED_DIR, svgfile)
            UPDATE_SVG && write(path, svg)
            @test svg == read(path, String)
        end
    end
end
