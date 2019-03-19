#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    DRAW_SETTING,
    BOND_DRAWER,
    atomnotation2d!,
    bondnotation2d!,
    boundary,
    chargesign,
    atommarkup,
    atomhtml


"""
    DRAW_SETTING

Default setting parameters of the molecule drawing canvas.

# Required fields

- `:display_terminal_carbon`(Bool) whether to display terminal C atom or not
- `:double_bond_notation`(Symbol)
    - `:alongside`: all double bonds are represented as a carbon skeleton and
                    a segment alongside it.
    - `:dual`: all double bonds are represented as two equal length parallel
               segments.
    - `:chain`: `:dual` for chain bonds and `:alongside` for ring bonds
    - `:terminal`: `:dual` for terminal bonds (adjacent to degree=1 node) and
                   `:alongside` for others (default)
- `:atomcolor`(Dict{Symbol,Color}) atom symbol and bond colors for organic atoms
- `:defaul_atom_color`(Dict{Symbol,Color}) colors for other atoms
"""
const DRAW_SETTING = Dict(
    :display_terminal_carbon => false,
    :double_bond_notation => :terminal,  # :alongside, :dual, :chain, :terminal
    :atomcolor => Dict(
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
    ),
    :default_atom_color => Color(0, 192, 192)
)


const BOND_DRAWER = Dict(
    1 => Dict(
        0 => singlebond!,
        1 => wedged!,
        6 => dashedwedged!,
        4 => wavesingle!
    ),
    2 => Dict(
        0 => clockwisedouble!,
        1 => counterdouble!,
        2 => doublebond!,
        3 => crossdouble!
    ),
    3 => Dict(
        0 => triplebond!
    )
)


abstract type Canvas end


"""
    atomnotation2d!(mol::VectorMol; setting=DRAW_SETTING)

Set atom symbol notation vectors of the molecule drawing
"""
function atomnotation2d!(mol::VectorMol; setting=DRAW_SETTING)
    haskey(mol, :AtomColor) && return
    elemental!(mol)
    atomc = setting[:atomcolor]
    dfc =  setting[:default_atom_color]
    mol[:Color2D] = Color[get(atomc, sym, dfc) for sym in mol[:Symbol]]
    termc = setting[:display_terminal_carbon]
    mol[:Visible2D] = Bool[
        (deg == 0 || sym != :C || (termc && deg == 1))
        for (deg, sym) in mol[:Degree, :Symbol]
    ]
    return
end


"""
    bondnotation2d!(mol::VectorMol; setting=DRAW_SETTING)

Set bond notation vectors of the molecule drawing
"""
function bondnotation2d!(mol::VectorMol; setting=DRAW_SETTING)
    haskey(mol, :BondNotation) && return
    haskey(mol.coords, :Cartesian2D) || throw(
        ErrorException(":Cartesian2D is required"))
    elemental!(mol)
    topology!(mol)
    if nodetype(mol) === SDFileAtom
        mol[:BondNotation] = getproperty.(edgevector(mol), :notation)
    else
        mol[:BondNotation] = zeros(Int, bondcount(mol))
    end
    dbnt = setting[:double_bond_notation]

    # All double bonds to be "="
    if dbnt == :dual
        for a in findall(mol[:BondOrder] .== 2)
            mol[:BondNotation][bond] = 2
        end
        return
    end

    # Or only non-ring bonds to be "="
    if dbnt == :terminal
        for a in findall(mol[:Degree] .== 1)
            bond = pop!(neighboredgeset(mol, a))
            if mol[:BondOrder][bond] == 2
                mol[:BondNotation][bond] = 2
            end
        end
    elseif dbnt == :chain
        for a in findall(.!mol[:RingBond] .* (mol[:BondOrder] .== 2))
            mol[:BondNotation][bond] = 2
        end
    end

    # Align double bonds alongside the ring
    rings = mol.annotation[:Topology].rings
    for ring in sort(rings, by=length, rev=true)
        cw = isclockwise(cyclicpath(mol.coords[:Cartesian2D], ring))
        if cw === nothing
            continue
        end
        succ = Dict()
        ordered = cw ? ring : reverse(ring)
        for n in 1:length(ring)-1
            succ[ordered[n]] = ordered[n+1]
        end
        succ[ordered[end]] = ordered[1]
        rsub = nodesubgraph(mol.graph, ring)
        for (i, e) in edgesiter(rsub)
            if mol[:BondOrder][i] == 2 && succ[e.u] != e.v
                mol[:BondNotation][i] = 1
            end
        end
    end
end


"""
    boundary(mol::VectorMol) -> (top, left, width, height, unit)

Get boundaries and an appropriate bond length unit for the molecule drawing
canvas.
"""
function boundary(mol::VectorMol)
    coords = mol.coords[:Cartesian2D]
    (left, right) = extrema(x_components(coords))
    (bottom, top) = extrema(y_components(coords))
    width = right - left
    height = top - bottom
    dists = []
    # Size unit
    for bond in edgevector(mol)
        u = _point(coords, bond.u)
        v =  _point(coords, bond.v)
        d = norm(v - u)
        if d > 0.0001  # Remove overlapped
            push!(dists, d)
        end
    end
    if isempty(dists)
        long = max(width, height)
        if long > 0.0001
            unit = long / sqrt(nodecount(mol))
        else
            unit = 1
        end
    else
        unit = median(dists) # Median bond length
    end

    return (top, left, width, height, unit)
end


"""
    chargesign(charge::Int) -> String

Get a charge sign.
"""
function chargesign(charge::Int)
    if charge == 0
        return ""
    end
    sign = charge > 0 ? "+" : "–" # en dash, not hyphen-minus
    num = abs(charge)
    return num > 1 ? string(num, sign) : sign
end


function atommarkup(symbol, charge, hcount, direction,
                    substart, subend, supstart, supend)
    if hcount == 1
        hs = "H"
    elseif hcount > 1
        hs = string("H", substart, hcount, subend)
    else
        hs = ""
    end
    chg = charge == 0 ? "" : string(supstart, chargesign(charge), supend)
    txt = direction == :left ? [chg, hs, symbol] : [symbol, hs, chg]
    return string(txt...)
end


function atomhtml(symbol, charge, hcount, direction)
    atommarkup(symbol, charge, hcount, direction,
               "<sub>", "</sub>", "<sup>", "</sup>")
end
