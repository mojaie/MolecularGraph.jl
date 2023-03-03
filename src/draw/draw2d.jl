#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    atom_color, is_atom_visible,
    single_bond_style, double_bond_style,
    chargesign, atommarkup, atomhtml,
    draw2d!, drawatomindex!, sethighlight!


"""
    atomcolor(mol::SimpleMolGraph; setting=DRAW_SETTING) -> Vector{Color}

Return atom colors for molecule 2D drawing
"""
atom_color(mol::SimpleMolGraph; color_theme=DEFAULT_ATOM_COLOR
    ) = [get(color_theme, sym, color_theme[:default]) for sym in atom_symbol(mol)]


"""
    isatomvisible(mol::SimpleMolGraph; setting=DRAW_SETTING) -> Vector{Bool}

Return whether the atom is visible in the 2D drawing.
"""
function is_atom_visible(mol::SimpleMolGraph; show_terminal_carbon=false)
    arr = (~).(init_node_descriptor(Bool, mol))
    deg_ = degree(mol)
    sym_ = atom_symbol(mol)
    chg_ = charge(mol)
    mul_ = multiplicity(mol)
    mas_ = getproperty.(vprops(mol), :mass)
    bondorder_ = bond_order(mol)
    for i in vertices(mol)
        sym_[i] === :C || continue
        chg_[i] == 0 || continue
        mul_[i] == 1 || continue
        mas_[i] === nothing || continue
        deg_[i] == 0 && continue
        deg_[i] == 1 && show_terminal_carbon && continue
        if deg_[i] == 2
            nbrs = neighbors(mol, i)
            u = edge_rank(mol, undirectededge(mol, i, nbrs[1]))
            v = edge_rank(mol, undirectededge(mol, i, nbrs[2]))
            if (bondorder_[u] == 2 && bondorder_[v] == 2)
                continue # allene-like
            end
        end
        arr[i] = false
    end
    return arr
end



function double_bond_style(mol::SimpleMolGraph)
    bondorder_ = bond_order(mol)
    coords = coords2d(mol)
    arr = init_edge_descriptor(Symbol, mol)
    for (i, e) in enumerate(edges(mol))
        if bondorder_[i] != 2
            arr[i] = :none
            continue
        end
        # TODO: other bond types which have notation
        if eproptype(mol) <: SDFBond && get_prop(mol, e, :notation) === 3
            arr[i] = :unspecified  # u x v (explicitly unspecified or racemic)
            continue
        end
        snbrs, dnbrs = edge_neighbors(mol, e)
        if length(snbrs) == 0 || length(dnbrs) == 0
            arr[i] = :none  # double bond at the end of chain
            continue
        end
        sdbs = map(snbrs) do snbr
            se = edge_rank(mol, undirectededge(mol, src(e), snbr))
            bondorder(mol)[se] == 2
        end
        ddbs = map(dnbrs) do dnbr
            de = edge_rank(mol, undirectededge(mol, dst(e), dnbr))
            bondorder(mol)[de] == 2
        end
        if any(sdbs) || any(ddbs)
            arr[i] = :none  # allene-like
            continue
        end
        arr[i] = :clockwise
    end
    # Align double bonds alongside the ring
    for ring in sort(sssr(mol), by=length, rev=true)
        cw = isclockwise(toarray(coords, ring))
        cw === nothing && continue
        ordered = cw ? ring : reverse(ring)
        rr = vcat(ordered, ordered)
        for i in 1:length(ordered)
            e = edge_rank(mol, undirectededge(mol, rr[i], rr[i + 1]))
            bondorder_[e] == 2 || continue
            arr[e] = rr[i] < rr[i + 1] ? :anticlockwise : :clockwise
        end
    end
    return arr
end


function single_bond_style(mol::SimpleMolGraph)
    # can be precalculated by coordgen
    has_prop(mol, :e_single_bond_style) && return get_prop(mol, :e_single_bond_style)
    bondorder = bond_order(mol)
    bondnotation = getproperty.(eprops(mol), :notation)
    isordered = getproperty.(eprops(mol), :isordered)
    arr = init_edge_descriptor(Symbol, mol)
    for (i, e) in enumerate(edges(mol))
        if bondorder[i] != 1
            arr[i] = :none
        elseif bondnotation[i] == 1
            arr[i] = isordered[i] ? :up : :revup
        elseif bondnotation[i] == 6
            arr[i] = isordered[i] ? :down : :revdown
        elseif bondnotation[i] == 4
            arr[i] = :unspecified
        else
            arr[i] = :none
        end
    end
    return arr
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
function draw2d!(canvas::Canvas, mol::SimpleMolGraph; kwargs...)
    # Canvas settings
    crds = coords2d(mol)
    initcanvas!(canvas, crds, boundary(mol, crds))
    canvas.valid || return

    # Properties
    atomsymbol_ = atom_symbol(mol)
    charge_ = charge(mol)
    implicith_ = implicit_hydrogens(mol)
    bondorder_ = bond_order(mol)
    atomcolor_ = atom_color(mol; kwargs...)
    isatomvisible_ = is_atom_visible(mol; kwargs...)
    bondstyle_ = map(
                single_bond_style(mol),
                double_bond_style(mol),
                bondorder_
            ) do sb, db, o
        if o == 1
            sb
        elseif o == 2
            db
        else
            :none
        end
    end

    # Draw bonds
    for (i, e) in enumerate(edges(mol))
        s, d = Tuple(e)
        setbond!(
            canvas, bondorder_[i], bondstyle_[i], s, d,
            atomcolor_[s], atomcolor_[d],
            isatomvisible_[s], isatomvisible_[d]
        )
    end

    # Draw atoms
    for i in vertices(mol)
        isatomvisible_[i] || continue
        pos = Point2D(canvas.coords, i)
        # Determine text direction
        if implicith_[i] > 0
            cosnbrs = []
            hrzn = pos + (1.0, 0.0)
            for nbr in neighbors(mol, i)
                posnbr = Point2D(canvas.coords, nbr)
                dist = distance(pos, posnbr)
                if dist > 0
                    dp = dot(hrzn - pos, posnbr - pos)
                    push!(cosnbrs, dp / dist)
                end
            end
            if isempty(cosnbrs) || minimum(cosnbrs) > 0
                # [atom]< or isolated node(ex. H2O, HCl)
                setatomright!(
                    canvas, pos, atomsymbol_[i], atomcolor_[i],
                    implicith_[i], charge_[i]
                )
                continue
            elseif maximum(cosnbrs) < 0
                # >[atom]
                setatomleft!(
                    canvas, pos, atomsymbol_[i], atomcolor_[i],
                    implicith_[i], charge_[i]
                )
                continue
            end
        end
        # -[atom]- or no hydrogens
        setatomcenter!(
            canvas, pos, atomsymbol_[i], atomcolor_[i],
            implicith_[i], charge_[i]
        )
    end
    return
end


function drawatomindex!(canvas::Canvas, mol::SimpleMolGraph;
                        color=Color(0, 0, 0), bgcolor=Color(240, 240, 255))
    isatomvisible_ = is_atom_visible(mol)
    for i in vertices(mol)
        offset = isatomvisible_[i] ? (0.0, canvas.fontsize/2.0) : (0.0, 0.0)
        pos = Point2D(canvas.coords, i) + offset
        setatomnote!(canvas, pos, string(i), color, bgcolor)
    end
    return
end


function sethighlight!(
        canvas::Canvas, substr::SimpleMolGraph; color=Color(253, 216, 53))
    isatomvisible_ = isatomvisible(substr)
    for e in edges(substr)
        setbondhighlight!(canvas, src(e), dst(e), color)
    end
    for i in vertices(substr)
        isatomvisible_[i] || continue
        pos = Point2D(canvas.coords, i)
        setatomhighlight!(canvas, pos, color)
    end
    return
end
