#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    DRAW_SETTING,
    atomcolor, isatomvisible, coords2d, bondnotation,
    chargesign, atommarkup, atomhtml,
    draw2d!, drawatomindex!


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


"""
    atomcolor(mol::VectorMol; setting=DRAW_SETTING)

Return atom colors for molecule 2D drawing
"""
function atomcolor(mol::VectorMol; setting=DRAW_SETTING)
    atomc = setting[:atomcolor]
    dfc =  setting[:default_atom_color]
    return [get(atomc, sym, dfc) for sym in mol[:atomsymbol]]
end


function isatomvisible(mol::VectorMol; setting=DRAW_SETTING)
    termc = setting[:display_terminal_carbon]
    vec = Bool[]
    for (deg, sym) in zip(mol[:nodedegree], mol[:atomsymbol])
        isvisible = deg == 0 || sym != :C || (termc && deg == 1)
        push!(vec, isvisible)
    end
    return vec
end


@cache function coords2d(mol::VectorMol; recalculate=false)
    recalculate && return compute2dcoords(mol)
    nodetype(mol) === SmilesAtom && return compute2dcoords(mol)
    matrix = zeros(Float64, nodecount(mol), 2)
    for (i, node) in nodesiter(mol)
        matrix[i, :] = node.coords[1:2]
    end
    return cartesian2d(matrix)
end


function bondnotation(mol::VectorMol; setting=DRAW_SETTING)
    if nodetype(mol) === SDFileAtom
        notation = getproperty.(edgevalues(mol), :notation)
    else
        notation = zeros(Int, edgecount(mol))
    end
    dbnt = setting[:double_bond_notation]

    # All double bonds to be "="
    if dbnt == :dual
        for b in findall(mol[:bondorder] .== 2)
            notation[b] = 2
        end
        return
    end

    # Or only non-ring bonds to be "="
    if dbnt == :terminal
        for i in findall(mol[:nodedegree] .== 1)
            b = pop!(incidences(mol, i))
            if mol[:bondorder][b] == 2
                notation[b] = 2
            end
        end
    elseif dbnt == :chain
        for b in findall(.!mol[:bond_isringmem] .* (mol[:bondorder] .== 2))
            notation[b] = 2
        end
    end

    # Align double bonds alongside the ring
    for ring in sort(sssr(mol), by=length, rev=true)
        cw = isclockwise(cyclicpath(mol[:coords2d], ring))
        cw === nothing && continue
        succ = Dict()
        ordered = cw ? ring : reverse(ring)
        for n in 1:length(ring)-1
            succ[ordered[n]] = ordered[n+1]
        end
        succ[ordered[end]] = ordered[1]
        rsub = nodesubgraph(mol.graph, ring)
        for (i, e) in edgesiter(rsub)
            if mol[:bondorder][i] == 2 && succ[e.u] != e.v
                notation[i] = 1
            end
        end
    end
    return notation
end


"""
    chargesign(charge::Int) -> String

Get a charge sign.
"""
function chargesign(charge::Int)
    charge == 0 && return ""
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


atomhtml(
    symbol, charge, hcount, direction
) = atommarkup(
    symbol, charge, hcount, direction, "<sub>", "</sub>", "<sup>", "</sup>")


"""
    draw2d!(canvas::Canvas, mol::VectorMol;
            setting=copy(DRAW_SETTING), recalculate=false)

Draw molecular image to the canvas.
"""
function draw2d!(canvas::Canvas, mol::VectorMol;
                 setting=copy(DRAW_SETTING), recalculate=false)
    # Canvas settings
    initcanvas!(canvas, mol)
    canvas.valid || return

    # Properties
    atomsymbol = mol[:atomsymbol]
    atomcolor = mol[:atomcolor]
    charge = mol[:charge]
    isatomvisible = mol[:isatomvisible]
    bondorder = mol[:bondorder]
    bondnotation = mol[:bondnotation]
    implicithcount = mol[:implicithcount]

    # Draw bonds
    for (i, bond) in edgesiter(mol)
        setbond!(
            canvas, bondorder[i], bondnotation[i],
            segment(canvas.coords, bond.u, bond.v),
            atomcolor[bond.u], atomcolor[bond.v],
            isatomvisible[bond.u], isatomvisible[bond.v]
        )
    end

    # Draw atoms
    for i in 1:atomcount(mol)
        isatomvisible[i] || continue
        pos = _point(canvas.coords, i)
        # Determine text direction
        if implicithcount[i] > 0
            cosnbrs = []
            hrzn = [pos[1] + 1.0, pos[2]]
            for nbr in adjacencies(mol, i)
                posnbr = _point(canvas.coords, nbr)
                dist = norm(posnbr - pos)
                if dist > 0
                    dp = dot(hrzn - pos, posnbr - pos)
                    push!(cosnbrs, dp / dist)
                end
            end
            if isempty(cosnbrs) || minimum(cosnbrs) > 0
                # [atom]< or isolated node(ex. H2O, HCl)
                setatomright!(
                    canvas, pos, atomsymbol[i], atomcolor[i],
                    implicithcount[i], charge[i]
                )
                continue
            elseif maximum(cosnbrs) < 0
                # >[atom]
                setatomleft!(
                    canvas, pos, atomsymbol[i], atomcolor[i],
                    implicithcount[i], charge[i]
                )
                continue
            end
        end
        # -[atom]- or no hydrogens
        setatomcenter!(
            canvas, pos, atomsymbol[i], atomcolor[i],
            implicithcount[i], charge[i]
        )
    end
    return
end


function drawatomindex!(canvas::Canvas, mol::VectorMol;
                        color=Color(0, 0, 0), bgcolor=Color(240, 240, 255))
    isatomvisible = mol[:isatomvisible]
    for i in 1:atomcount(mol)
        offset = isatomvisible[i] ? [0 canvas.fontsize/2] : [0 0]
        pos = _point(canvas.coords, i) + offset
        setatomnote!(canvas, pos, string(i), color, bgcolor)
    end
    return
end
