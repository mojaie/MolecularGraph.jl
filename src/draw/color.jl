#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# Default extended CPK theme for 2D atom symbol drawing

const DEFAULT_ATOM_COLOR = Dict(k => RGB{N0f8}((v ./ 255)...) for (k, v) in Dict(
    :default => (0, 192, 192),
    :H => (0, 0, 0),
    :B => (128, 0, 0),
    :C => (0, 0, 0),
    :N => (0, 0, 255),
    :O => (255, 0, 0),
    :F => (0, 255, 0),
    :Si => (128, 64, 192),
    :P => (192, 0, 192),
    :S => (192, 192, 0),
    :Cl => (64, 192, 64),
    :As => (128, 0, 128),
    :Se => (128, 128, 0),
    :Br => (0, 192, 0),
    :I => (0, 128, 0)
))


# https://www.umass.edu/microbio/rasmol/distrib/rasman.htm#cpkcolours
# default for 3D model

const RASMOL_ATOM_COLOR = Dict(k => RGB{N0f8}((v ./ 255)...) for (k, v) in Dict(
    :default => (255,20,147),
    :H => (255,255,255),
    :B => (0,255,0),
    :C => (200,200,200),
    :N => (143,143,255),
    :O => (240,0,0),
    :F => (218,165,32),
    :Si => (218,165,32),
    :P => (255,165,0),
    :S => (255,200,50),
    :Cl => (0,255,0),
    :Br => (165,42,42),
    :I => (160,32,240),
    :Ba => (255,165,0),
    :Fe => (255,165,0),
    :Na => (0,0,255),
    :Mg => (34,139,34),
    :Zn => (165,42,42),
    :Cu => (165,42,42),
    :Ni => (165,42,42),
    :Ca => (128,128,144),
    :Mn => (128,128,144),
    :Al => (128,128,144),
    :Ti => (128,128,144),
    :Cr => (128,128,144),
    :Ag => (128,128,144),
    :Au => (218,165,32),
    :Li => (178,34,34),
    :He => (255,192,203)
))


"""
    atom_color(mol::SimpleMolGraph) -> Vector{RGB}

Return atom colors for molecule 2D drawing
"""
atom_color(mol::SimpleMolGraph; kwargs...
    ) = [atom_color(mol[i]; kwargs...) for i in vertices(mol)]

atom_color(atom::StandardAtom; color_theme=DEFAULT_ATOM_COLOR, kwargs...
    ) = get(color_theme, atom_symbol(atom), color_theme[:default])

atom_color(atom::AbstractAtom; kwargs...) = atom_color(
    atom,
    Val(has_hydrogens(typeof(atom))),
    Val(has_label(typeof(atom)))
    ; kwargs...
)
atom_color(atom, ::Val{true}, ::Val{false}; kwargs...) = atom_color(atom.center)
atom_color(atom, ::Val, ::Val{true}; color_theme=DEFAULT_ATOM_COLOR, kwargs...) = color_theme[:C]


"""
    atom_coloralpha(mol::SimpleMolGraph; alpha=1.0) -> Vector{RGBA}

Return atom colors for molecule 2D drawing
"""
atom_coloralpha(mol::SimpleMolGraph; alpha=1.0, kwargs...
    ) = [coloralpha(c, alpha) for c in atom_color(mol; kwargs...)]
