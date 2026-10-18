# Drawing molecules

This page deals with 2D/3D rendering of molecules.

- 2D structure images
  - SVG images
  - Drawing options
  - PNG images
  - HTML layout
  - Regenerating 2D coordinates
- 3D molecule rendering using Makie.jl

2D SVG drawing is available out of the box. Some features depend on optional packages
that are loaded as package extensions:

- `Cairo` is required to export PNG images of 2D structures.
- `Makie` (with a backend such as `GLMakie` or `CairoMakie`) is required to render 3D
  molecular models. Note that the loading time of Makie can be substantially longer.

```@example drawing
using MolecularGraph, Graphs
```


## 2D structure images

### SVG images

In this page, we use a molecule loaded from an SDFile. The SDFiles used here can be
obtained from PubChem (e.g. `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/6437877/SDF`
for cefditoren pivoxil).

```@example drawing
datadir = joinpath(pkgdir(MolecularGraph), "docs", "src", "assets", "data")
mol = sdftomol(joinpath(datadir, "Cefditoren Pivoxil.mol"))
remove_hydrogens!(mol)
nothing # hide
```

[`drawsvg`](@ref) returns the SVG text of the molecular structure.

```@example drawing
svg = drawsvg(mol)
print(first(svg, 300), "...")
```

To display the SVG in an HTML environment (Documenter, Pluto, Jupyter...), wrap the string
in `HTML`. The SVG text can also be saved to a file directly
(e.g. `write("mol.svg", drawsvg(mol))`).

The default SVG width and height are "100%", so the actual size is usually determined by
the container HTML element. The size can also be specified by the `width` and `height` options.

```@example drawing
HTML(drawsvg(mol, width=200, height=200))
```

### Drawing options

- `atomindex=true` displays atom indices (`indexcolor` and `indexbgcolor` change the colors).
- `atomhighlight` (a vector of atom indices) and `bondhighlight` (a vector of `Edge`s) highlight
  components of interest. `bondhighlight` should contain only existing bonds of the molecule,
  otherwise an `ArgumentError` is thrown. The color can be changed by `highlightcolor`.

```@example drawing
HTML(drawsvg(mol, width=300, height=300, atomindex=true))
```

```@example drawing
HTML(drawsvg(mol, width=300, height=300,
    atomhighlight=[14, 15], bondhighlight=[Edge(4 => 19), Edge(7 => 25)]))
```

Colors are given as strings of hex codes (e.g. `"#FDD835"`, `"#FFF"`) or CSS color names
(e.g. `"lightblue"`). The background is set by `bgcolor` and `bgopacity`.

```@example drawing
HTML(drawsvg(mol, width=300, height=300,
    bgcolor="#F5F5F5", highlightcolor="lightblue",
    bondhighlight=Edge.([(3, 39), (14, 39), (14, 31), (31, 37), (3, 37)])))
```

By default, carbon atom labels are not shown. The `show_carbon` option changes this behavior:
`:simple` (default) hides carbon labels, `:terminal` shows only labels of terminal carbons,
and `:all` shows all carbon labels.

```@example drawing
HTML(drawsvg(smilestomol("CC(C)CC(=O)O"), width=200, height=200, show_carbon=:all))
```

### PNG images

[`drawpng`](@ref) writes PNG image data into an IO object with the given `width` and `height`
(in px) and drawing options similar to those of `drawsvg`. This requires `Cairo` to be loaded.

```julia
using MolecularGraph, Cairo

open("mol.png", "w") do io
    drawpng(io, mol, 300, 300; atomindex=true)
end
```

### HTML layout

[`html_fixed_size`](@ref) is a convenient function to wrap an SVG image in an HTML element
of a fixed size. The keyword arguments are passed to `drawsvg`.

```@example drawing
html_fixed_size(mol, 150, 150,
    bondhighlight=Edge.([(3, 39), (14, 39), (14, 31), (31, 37), (3, 37)]))
```

In HTML environments, a molecule object as the result of an expression is displayed by
`html_fixed_size(mol, "250px", "250px")` by default.

[`html_grid`](@ref)`(mols, cols, rowheight; kwargs...)` displays SVG images of molecules in a grid
with `cols` columns and the given row height. The width is determined by the parent HTML
element. A vector of molecules is displayed by `html_grid(mols, 4, "200px")` by default.

```@example drawing
mols = [
    smilestomol(raw"O=S(=O)(N)c1c(Cl)cc(c(C(=O)O)c1)NCc2occc2"),
    smilestomol(raw"CCC(CC)OC1C=C(CC(C1NC(=O)C)N)C(=O)OCC"),
    smilestomol(raw"O=C([C@H](CCC(O)=O)NC(C1=CC=C(N(CC2=CN=C(N=C(N)N=C3N)C3=N2)C)C=C1)=O)O"),
    smilestomol(raw"O=C2\N=C(/OC2c1ccccc1)N"),
    smilestomol(raw"O(CCN(C)C)C(c1ccccc1)c2ccccc2"),
    smilestomol(raw"C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O")
]
html_grid(mols, 3, "250px", atomindex=true)
```

### Regenerating 2D coordinates

As SMILES does not have atom coordinates, 2D coordinates are generated when the molecule is
drawn. Internally, MolecularGraph.jl uses Schrodinger's coordgenlibs
([https://github.com/schrodinger/coordgenlibs](https://github.com/schrodinger/coordgenlibs))
for 2D coordinate generation.

```@example drawing
smol = smilestomol(raw"O=C3N2/C(=C(/C=C\c1scnc1C)CS[C@@H]2[C@@H]3NC(=O)C(=N\OC)/c4nc(sc4)N)C(=O)O")
html_fixed_size(smol, 250, 250, atomindex=true)
```

[`coordgen!`](@ref) (re)generates 2D coordinates and stores them in the molecule, replacing
the existing ones. This also works with molecules from SDFiles, so you can recalculate 2D
coordinates after modifying the molecular graph, or replace the original coordinates of
the SDFile. [`coordgen`](@ref) returns the coordinates and bond styles without modifying the molecule.

```@example drawing
mol2 = sdftomol(joinpath(datadir, "Cefditoren Pivoxil.mol"))
remove_hydrogens!(mol2)
coordgen!(mol2)
html_grid([mol, mol2], 2, "250px")
```


## 3D molecule rendering using Makie.jl

3D molecular rendering is implemented with Makie.jl. This requires one of the Makie backends
(e.g. `GLMakie`) to be loaded. 3D coordinates are taken from the SDFile (e.g. PubChem 3D
conformers, `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/54670067/SDF?record_type=3d`).

Available rendering functions are [`spacefilling`](@ref), [`ballstick`](@ref), [`stick`](@ref)
and [`wire`](@ref), and their in-place versions `spacefilling!`, `ballstick!`, `stick!` and
`wire!`, which can be used similarly to Makie's `meshscatter!` and `mesh!`.

```julia
using MolecularGraph, GLMakie

mol3d = sdftomol(joinpath(datadir, "ascorbic_acid_3d.mol"))

fig = Figure(figure_padding=0)
ax = LScene(fig[1, 1], show_axis=false, width=300, height=300)
spacefilling!(mol3d)
resize_to_layout!(fig)
fig
```

Main keyword options:

- `spacefilling`: `radii` (default `"van der Waals"`)
- `ballstick`: `radii` (ball size; a number, or a `Dict` of atomic number to radius) and `bonddiameter`
- `stick` and `wire`: `size` (diameter of the sticks/wires)

```julia
fig = Figure(figure_padding=0)
ax = LScene(fig[1, 1], show_axis=false, width=300, height=300)
plt = ballstick!(mol3d, radii=Dict(1 => 0.25, 6 => 0.4, 8 => 0.4))
resize_to_layout!(fig)
fig
```

The returned plot object can be manipulated with Makie functions (e.g. `rotate!`) to create
animations with `record`.
