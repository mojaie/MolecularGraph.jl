#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    DRAW_SETTING,
    atomcolor, isatomvisible, sdfcoords2d, coords2d,
    chargesign, atommarkup, atomhtml,
    draw2d!, drawatomindex!, sethighlight!


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
    :double_bond_style => :terminal,  # :alongside, :dual, :chain, :terminal
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

# For 3d rendering, we use white for hydrogen
const DRAW_SETTING3 = (dct = deepcopy(DRAW_SETTING); dct[:atomcolor][:H] = Color(255, 255, 255); dct)

"""
    atomcolor(mol::GraphMol; setting=DRAW_SETTING) -> Vector{Color}

Return atom colors for molecule 2D drawing
"""
@cachefirst function atomcolor(
        mol::GraphMol; kwargs...)
    atomcolor(atomsymbol(mol); kwargs...)
end

function atomcolor(syms::AbstractVector; setting=DRAW_SETTING, kwargs...)
    atomc = setting[:atomcolor]
    dfc = setting[:default_atom_color]
    return [get(atomc, sym, dfc) for sym in syms]
end

atomcolor(view::SubgraphView; kwargs...) = atomcolor(view.graph; kwargs...)


"""
    isatomvisible(mol::GraphMol; setting=DRAW_SETTING) -> Vector{Bool}

Return whether the atom is visible in the 2D drawing.
"""
@cachefirst function isatomvisible(
        mol::GraphMol; setting=DRAW_SETTING, kwargs...)
    termc = setting[:display_terminal_carbon]
    vec = Bool[]
    deg_ = nodedegree(mol)
    sym_ = atomsymbol(mol)
    for i in 1:nodecount(mol)
        isvisible = deg_[i] == 0 || sym_[i] !== :C || (termc && deg_[i] == 1)
        push!(vec, isvisible)
    end
    return vec
end

isatomvisible(view::SubgraphView; kwargs...
    ) = isatomvisible(view.graph; kwargs...)



@cachefirst function sdfcoords2d(mol::SDFile)
    coords = zeros(Float64, nodecount(mol), 2)
    for (i, node) in enumerate(nodeattrs(mol))
        coords[i, :] = node.coords[1:2]
    end
    return coords
end

sdfcoords2d(view::SubgraphView) = sdfcoords2d(view.graph)



@cachefirst function coords2d(
        mol::GraphMol; forcecoordgen=false, setting=DRAW_SETTING, kwargs...)
    if nodeattrtype(mol) === SDFileAtom && !forcecoordgen
        coords = sdfcoords2d(mol)
        style = getproperty.(edgeattrs(mol), :notation)
    else
        coords, style = coordgen(mol)
    end
    bondorder_ = bondorder(mol)

    # All double bonds to be "="
    if setting[:double_bond_style] == :dual
        for b in findall(bondorder_ .== 2)
            style[b] = 2
        end
        return
    end

    # Or only non-ring bonds to be "="
    if setting[:double_bond_style] == :terminal
        for i in findall(nodedegree(mol) .== 1)
            b = iterate(incidences(mol, i))[1]
            if bondorder_[b] == 2
                style[b] = 2
            end
        end
    elseif setting[:double_bond_style] == :chain
        for b in findall(.!isringbond(mol) .* (bondorder_ .== 2))
            style[b] = 2
        end
    end

    # Align double bonds alongside the ring
    for ring in sort(sssr(mol), by=length, rev=true)
        cw = isclockwise(toarray(coords, ring))
        cw === nothing && continue
        ordered = cw ? ring : reverse(ring)
        rr = vcat(ordered, ordered)
        for i in 1:length(ordered)
            e = findedgekey(mol, rr[i], rr[i + 1])
            (u, v) = getedge(mol, e)
            if bondorder_[e] == 2 && u != rr[i]
                style[e] = 1
            end
        end
    end

    return coords, style
end

coords2d(view::SubgraphView; kwargs...) = coords2d(view.graph; kwargs...)



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


function atommarkup(symbol, charge, hydrogenconnected, direction,
                    substart, subend, supstart, supend)
    if hydrogenconnected == 1
        hs = "H"
    elseif hydrogenconnected > 1
        hs = string("H", substart, hydrogenconnected, subend)
    else
        hs = ""
    end
    chg = charge == 0 ? "" : string(supstart, chargesign(charge), supend)
    txt = direction == :left ? [chg, hs, symbol] : [symbol, hs, chg]
    return string(txt...)
end


atomhtml(
    symbol, charge, hydrogenconnected, direction
) = atommarkup(
    symbol, charge, hydrogenconnected, direction, "<sub>", "</sub>", "<sup>", "</sup>")


"""
    draw2d!(canvas::Canvas, mol::UndirectedGraph; kwargs...)

Draw molecular image to the canvas.
"""
function draw2d!(canvas::Canvas, mol::UndirectedGraph; kwargs...)
    # Canvas settings
    coords_, bondstyles_ = coords2d(mol; kwargs...)
    initcanvas!(canvas, coords_, boundary(mol, coords_))
    canvas.valid || return

    # Properties
    atomsymbol_ = atomsymbol(mol)
    atomcolor_ = atomcolor(mol; kwargs...)
    charge_ = charge(mol)
    isatomvisible_ = isatomvisible(mol; kwargs...)
    bondorder_ = bondorder(mol)
    implicithconnected_ = implicithconnected(mol)

    # Draw bonds
    for i in edgeset(mol)
        (u, v) = getedge(mol, i)
        if bondorder_[i] == 1 && bondstyles_[i] in (2, 7)
            # cordgen reverse bonds
            (t, s) = (u, v)
            bs = bondstyles_[i] - 1
        else
            (s, t) = (u, v)
            bs = bondstyles_[i]
        end
        setbond!(
            canvas, bondorder_[i], bs,
            Segment{Point2D}(canvas.coords, s, t),
            atomcolor_[u], atomcolor_[v],
            isatomvisible_[u], isatomvisible_[v]
        )
    end

    # Draw atoms
    for i in nodeset(mol)
        isatomvisible_[i] || continue
        pos = Point2D(canvas.coords, i)
        # Determine text direction
        if implicithconnected_[i] > 0
            cosnbrs = []
            hrzn = pos + (1.0, 0.0)
            for adj in adjacencies(mol, i)
                posnbr = Point2D(canvas.coords, adj)
                dist = Geometry.distance(pos, posnbr)
                if dist > 0
                    dp = dot(hrzn - pos, posnbr - pos)
                    push!(cosnbrs, dp / dist)
                end
            end
            if isempty(cosnbrs) || minimum(cosnbrs) > 0
                # [atom]< or isolated node(ex. H2O, HCl)
                setatomright!(
                    canvas, pos, atomsymbol_[i], atomcolor_[i],
                    implicithconnected_[i], charge_[i]
                )
                continue
            elseif maximum(cosnbrs) < 0
                # >[atom]
                setatomleft!(
                    canvas, pos, atomsymbol_[i], atomcolor_[i],
                    implicithconnected_[i], charge_[i]
                )
                continue
            end
        end
        # -[atom]- or no hydrogens
        setatomcenter!(
            canvas, pos, atomsymbol_[i], atomcolor_[i],
            implicithconnected_[i], charge_[i]
        )
    end
    return
end


function drawatomindex!(canvas::Canvas, mol::UndirectedGraph;
                        color=Color(0, 0, 0), bgcolor=Color(240, 240, 255))
    isatomvisible_ = isatomvisible(mol)
    for i in nodeset(mol)
        offset = isatomvisible_[i] ? (0.0, canvas.fontsize/2.0) : (0.0, 0.0)
        pos = Point2D(canvas.coords, i) + offset
        setatomnote!(canvas, pos, string(i), color, bgcolor)
    end
    return
end


function sethighlight!(
        canvas::Canvas, substr::UndirectedGraph; color=Color(253, 216, 53))
    isatomvisible_ = isatomvisible(substr)
    for i in edgeset(substr)
        (u, v) = getedge(substr, i)
        setbondhighlight!(canvas, Segment{Point2D}(canvas.coords, u, v), color)
    end
    for i in nodeset(substr)
        isatomvisible_[i] || continue
        pos = Point2D(canvas.coords, i)
        setatomhighlight!(canvas, pos, color)
    end
    return
end
