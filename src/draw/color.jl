#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# TODO: common draw color theme based on Colors.jl

struct Color
    r::Int
    g::Int
    b::Int
end

# Default extended CPK theme for 2D atom symbol drawing

const DEFAULT_ATOM_COLOR = Dict(
    :default => Color(0, 192, 192),
    :H => Color(0, 0, 0),
    :B => Color(128, 0, 0),
    :C => Color(0, 0, 0),
    :N => Color(0, 0, 255),
    :O => Color(255, 0, 0),
    :F => Color(0, 255, 0),
    :Si => Color(128, 64, 192),
    :P => Color(192, 0, 192),
    :S => Color(192, 192, 0),
    :Cl => Color(64, 192, 64),
    :As => Color(128, 0, 128),
    :Se => Color(128, 128, 0),
    :Br => Color(0, 192, 0),
    :I => Color(0, 128, 0)
)


# https://www.umass.edu/microbio/rasmol/distrib/rasman.htm#cpkcolours
# default for 3D model

const RASMOL_ATOM_COLOR = Dict(
    :default => Color(255,20,147),
    :H => Color(255,255,255),
    :B => Color(0,255,0),
    :C => Color(200,200,200),
    :N => Color(143,143,255),
    :O => Color(240,0,0),
    :F => Color(218,165,32),
    :Si => Color(218,165,32),
    :P => Color(255,165,0),
    :S => Color(255,200,50),
    :Cl => Color(0,255,0),
    :Br => Color(165,42,42),
    :I => Color(160,32,240),
    :Ba => Color(255,165,0),
    :Fe => Color(255,165,0),
    :Na => Color(0,0,255),
    :Mg => Color(34,139,34),
    :Zn => Color(165,42,42),
    :Cu => Color(165,42,42),
    :Ni => Color(165,42,42),
    :Ca => Color(128,128,144),
    :Mn => Color(128,128,144),
    :Al => Color(128,128,144),
    :Ti => Color(128,128,144),
    :Cr => Color(128,128,144),
    :Ag => Color(128,128,144),
    :Au => Color(218,165,32),
    :Li => Color(178,34,34),
    :He => Color(255,192,203)
)


"""
    atomcolor(mol::SimpleMolGraph; setting=DRAW_SETTING) -> Vector{Color}

Return atom colors for molecule 2D drawing
"""
atom_color(mol::SimpleMolGraph; color_theme=DEFAULT_ATOM_COLOR, kwargs...
    ) = [get(color_theme, sym, color_theme[:default]) for sym in atom_symbol(mol)]