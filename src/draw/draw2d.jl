#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


"""
    is_atom_visible(mol::SimpleMolGraph; setting=DRAW_SETTING) -> Vector{Bool}

Return whether the atom is visible in the 2D drawing.
`show_carbon`: `simple` - does not show carbon labels, `terminal` - show only terminal carbonlabels,
`all` - show all carbon labels.
"""
function is_atom_visible(g, sym, chg, mul, ms, bo; show_carbon=:simple, kwargs...)
    arr = trues(nv(g))
    show_carbon === :all && return arr
    er = Dict(e => i for (i, e) in enumerate(edges(g)))
    for i in vertices(g)
        sym[i] === :C || continue
        chg[i] == 0 || continue
        mul[i] == 1 || continue
        ms[i] === nothing || continue
        degree(g, i) == 0 && continue
        degree(g, i) == 1 && show_carbon === :terminal && continue
        if degree(g, i) == 2
            nbrs = neighbors(g, i)
            u = er[undirectededge(g, i, nbrs[1])]
            v = er[undirectededge(g, i, nbrs[2])]
            if (bo[u] == 2 && bo[v] == 2)
                continue # allene-like
            end
        end
        arr[i] = false
    end
    return arr
end

function is_atom_visible(mol::SimpleMolGraph; show_carbon=:simple, kwargs...)
    get_state(mol, :has_updates) && dispatch!(mol, :on_update)
    return is_atom_visible(mol.graph, atom_symbol(mol), charge(mol), multiplicity(mol),
        [get_prop(mol, i, :mass) for i in vertices(mol)], bond_order(mol); kwargs...)
end


function double_bond_style(g, bondorder_, ntt, coords, sssr_)
    arr = Vector{Symbol}(undef, ne(g))
    er = Dict(e => i for (i, e) in enumerate(edges(g)))
    for (i, e) in enumerate(edges(g))
        if bondorder_[i] != 2
            arr[i] = :none
            continue
        end
        if ntt[i] == 3
            arr[i] = :unspecified  # u x v (explicitly unspecified or racemic)
            continue
        end
        snbrs, dnbrs = edge_neighbors(g, e)
        if length(snbrs) == 0 || length(dnbrs) == 0
            arr[i] = :none  # double bond at the end of chain
            continue
        end
        sdbs = map(snbrs) do snbr
            se = er[undirectededge(g, src(e), snbr)]
            bondorder_[se] == 2
        end
        ddbs = map(dnbrs) do dnbr
            de = er[undirectededge(g, dst(e), dnbr)]
            bondorder_[de] == 2
        end
        if any(sdbs) || any(ddbs)
            arr[i] = :none  # allene-like
            continue
        end
        arr[i] = :clockwise
    end
    # Align double bonds alongside the ring
    for ring in sort(sssr_, by=length, rev=true)
        cw = isclockwise(toarray(coords, ring))
        cw === nothing && continue
        ordered = cw ? ring : reverse(ring)
        rr = vcat(ordered, ordered)
        for i in 1:length(ordered)
            e = er[undirectededge(g, rr[i], rr[i + 1])]
            bondorder_[e] == 2 || continue
            arr[e] = rr[i] < rr[i + 1] ? :clockwise : :anticlockwise
        end
    end
    return arr
end

function double_bond_style(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :on_update)
    ntt = hasfield(eproptype(mol), :notation
        ) ? [get_prop(mol, e, :notation) for e in edges(mol)] : zeros(ne(mol))
    return double_bond_style(mol.graph, bond_order(mol), ntt, coords2d(mol), sssr(mol))
end


function single_bond_style(bondorder, bondnotation, isordered)
    arr = Vector{Symbol}(undef, length(bondorder))
    for i in 1:length(bondorder)
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

function single_bond_style(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :on_update)
    # can be precalculated by coordgen
    has_state(mol, :e_single_bond_style) && return get_state(mol, :e_single_bond_style)
    return single_bond_style(
        bond_order(mol), [get_prop(mol, e, :notation) for e in edges(mol)],
        [get_prop(mol, e, :isordered) for e in edges(mol)]
    )
end


function bond_style(bondorder, singlebondstyle, doublebondstyle)
    arr = Vector{Symbol}(undef, length(bondorder))
    for i in 1:length(bondorder)
        if bondorder[i] == 1
            arr[i] = singlebondstyle[i]
        elseif bondorder[i] == 2
            arr[i] = doublebondstyle[i]
        else
            arr[i] = :none
        end
    end
    return arr
end

bond_style(mol::SimpleMolGraph) = bond_style(bond_order(mol), single_bond_style(mol), double_bond_style(mol))



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
    # get coords
    if !hasfield(vproptype(mol), :coords) && !has_state(mol, :v_coords2d)  # default SMILESAtom
        crds, sb_style = coordgen(mol)
    else  # SDFAtom or has coordgen! precache
        crds = coords2d(mol)
        sb_style = single_bond_style(mol)
    end
    # Canvas settings
    initcanvas!(canvas, crds, boundary(mol, crds))
    canvas.valid || return

    # Properties
    atomsymbol_ = atom_symbol(mol)
    charge_ = charge(mol)
    implicith_ = implicit_hydrogens(mol)
    bondorder_ = bond_order(mol)
    atomcolor_ = atom_color(mol; kwargs...)
    isatomvisible_ = is_atom_visible(mol; kwargs...)
    bondstyle_ = bond_style(bondorder_, sb_style, double_bond_style(mol))

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


function drawatomindex!(canvas::Canvas, isatomvisible, color, bgcolor)
    for (i, v) in enumerate(isatomvisible)
        offset = v ? (0.0, canvas.fontsize/2.0) : (0.0, 0.0)
        pos = Point2D(canvas.coords, i) + offset
        setatomnote!(canvas, pos, string(i), color, bgcolor)
    end
end


function sethighlight!(canvas::Canvas, edge_list::Vector{<:Edge}, color)
    for e in edge_list
        setbondhighlight!(canvas, src(e), dst(e), color)
    end
end

function sethighlight!(canvas::Canvas, node_list::Vector{<:Integer}, color)
    for i in node_list
        pos = Point2D(canvas.coords, i)
        setatomhighlight!(canvas, pos, color)
    end
end
