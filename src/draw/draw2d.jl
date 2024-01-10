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
            u = er[u_edge(g, i, nbrs[1])]
            v = er[u_edge(g, i, nbrs[2])]
            if (bo[u] == 2 && bo[v] == 2)
                continue # allene-like
            end
        end
        arr[i] = false
    end
    return arr
end

function is_atom_visible(mol::SimpleMolGraph; show_carbon=:simple, kwargs...)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    return is_atom_visible(mol.graph, atom_symbol(mol), charge(mol), multiplicity(mol),
        [get_prop(mol, i, :mass) for i in vertices(mol)], bond_order(mol); kwargs...)
end


function double_bond_style(g, bondorder_, coords, sssr_)
    arr = Vector{Symbol}(undef, ne(g))
    er = Dict(e => i for (i, e) in enumerate(edges(g)))
    for (i, e) in enumerate(edges(g))
        if bondorder_[i] != 2
            arr[i] = :none
            continue
        end
        snbrs, dnbrs = edge_neighbors(g, e)
        if length(snbrs) == 0 || length(dnbrs) == 0
            arr[i] = :none  # double bond at the end of chain
            continue
        end
        sdbs = map(snbrs) do snbr
            se = er[u_edge(g, src(e), snbr)]
            bondorder_[se] == 2
        end
        ddbs = map(dnbrs) do dnbr
            de = er[u_edge(g, dst(e), dnbr)]
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
            e = er[u_edge(g, rr[i], rr[i + 1])]
            bondorder_[e] == 2 || continue
            arr[e] = rr[i] < rr[i + 1] ? :clockwise : :anticlockwise
        end
    end
    return arr
end


function sdf_bond_style(bondorder, bondnotation, isordered)
    arr = Vector{Symbol}(undef, length(bondorder))
    for i in 1:length(bondorder)
        if bondnotation[i] == 3
            arr[i] = :cis_trans
        elseif bondorder[i] != 1
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

function sdf_bond_style(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    return sdf_bond_style(
        bond_order(mol),
        [get_prop(mol, e, :notation) for e in edges(mol)],
        [get_prop(mol, e, :isordered) for e in edges(mol)]
    )
end


function bond_style(g, bondorder, defaultbondstyle, coords_, sssr_)
    doublebondstyle = double_bond_style(g, bondorder, coords_, sssr_)
    arr = copy(defaultbondstyle)
    for i in 1:length(bondorder)
        defaultbondstyle[i] === :cis_trans && continue
        if bondorder[i] == 2
            arr[i] = doublebondstyle[i]
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


function atommarkup(symbol, charge, implicith, direction,
                    substart, subend, supstart, supend)
    if implicith == 1
        hs = "H"
    elseif implicith > 1
        hs = string("H", substart, implicith, subend)
    else
        hs = ""
    end
    chg = charge == 0 ? "" : string(supstart, chargesign(charge), supend)
    txt = direction == :left ? [chg, hs, symbol] : [symbol, hs, chg]
    return string(txt...)
end


atomhtml(
    symbol, charge, implicith, direction
) = atommarkup(
    symbol, charge, implicith, direction, "<sub>", "</sub>", "<sup>", "</sup>")


"""
    boundary(mol::SimpleMolGraph, coords::AbstractArray{Float64}
        ) -> (top, left, width, height, unit)

Get boundaries and an appropriate bond length unit for the molecule drawing
canvas.
"""
function boundary(mol::SimpleMolGraph, coords::AbstractArray{Float64})
    (left, right) = extrema(x_components(coords))
    (bottom, top) = extrema(y_components(coords))
    width = right - left
    height = top - bottom
    dists = []
    # Size unit
    for e in edges(mol)
        d = distance(Point2D(coords, src(e)), Point2D(coords, dst(e)))
        if d > 0.0001  # Remove overlapped
            push!(dists, d)
        end
    end
    if isempty(dists)
        long = max(width, height)
        unit = long > 0.0001 ? long / sqrt(nv(mol)) : 1
    else
        unit = median(dists) # Median bond length
    end
    return (top, left, width, height, unit)
end


function trimbond(canvas::Canvas, seg, uvis, vvis)
    if uvis && vvis
        return trim_uv(seg, canvas.trimoverlapf * 2)
    elseif uvis
        return trim_u(seg, canvas.trimoverlapf)
    elseif vvis
        return trim_v(seg, canvas.trimoverlapf)
    end
    return seg
end

function singlebond!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis)
    seg = trimbond(canvas, Segment{Point2D}(canvas.coords, u, v), uvis, vvis)
    drawline!(canvas, seg, ucolor, vcolor)
    return
end

function wedged!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis)
    seg = trimbond(canvas, Segment{Point2D}(canvas.coords, u, v), uvis, vvis)
    drawwedge!(canvas, seg, ucolor, vcolor)
    return
end

function reversed_wedged!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis)
    seg = trimbond(canvas, Segment{Point2D}(canvas.coords, v, u), vvis, uvis)
    drawwedge!(canvas, seg, vcolor, ucolor)
    return
end

function dashedwedged!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis)
    seg = trimbond(canvas, Segment{Point2D}(canvas.coords, u, v), uvis, vvis)
    drawdashedwedge!(canvas, seg, ucolor, vcolor)
    return
end

function reversed_dashedwedged!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis)
    seg = trimbond(canvas, Segment{Point2D}(canvas.coords, v, u), vvis, uvis)
    drawdashedwedge!(canvas, seg, vcolor, ucolor)
    return
end

function wavesingle!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis)
    seg = trimbond(canvas, Segment{Point2D}(canvas.coords, u, v), uvis, vvis)
    drawwave!(canvas, seg, ucolor, vcolor)
    return
end

function doublebond!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis)
    dist = canvas.scaleunit * canvas.mbwidthf / 2
    seg = trimbond(canvas, Segment{Point2D}(canvas.coords, u, v), uvis, vvis)
    seg1 = translate(seg, pi / 2, dist)
    seg2 = translate(seg, -pi / 2, dist)
    drawline!(canvas, seg1, ucolor, vcolor)
    drawline!(canvas, seg2, ucolor, vcolor)
    return
end

function crossdouble!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis)
    dist = canvas.scaleunit * canvas.mbwidthf / 2
    seg = trimbond(canvas, Segment{Point2D}(canvas.coords, u, v), uvis, vvis)
    cw = translate(seg, pi / 2, dist)
    ccw = translate(seg, -pi / 2, dist)
    drawline!(canvas, Segment(cw.u, ccw.v), ucolor, vcolor)
    drawline!(canvas, Segment(ccw.u, cw.v), ucolor, vcolor)
    return
end

function ringdouble!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis, direction)
    seg = trimbond(canvas, Segment{Point2D}(canvas.coords, u, v), uvis, vvis)
    dist = canvas.scaleunit * canvas.mbwidthf
    segin = trim_uv(translate(seg, direction, dist), canvas.triminnerf)
    drawline!(canvas, seg, ucolor, vcolor)
    drawline!(canvas, segin, ucolor, vcolor)
    return
end

# NOTE: the direction is reversed due to x-axis reflection
clockwisedouble!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis
    ) = ringdouble!(canvas, u, v, ucolor, vcolor, uvis, vvis, pi / 2)
counterdouble!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis
    ) = ringdouble!(canvas, u, v, ucolor, vcolor, uvis, vvis, -pi / 2)

function triplebond!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis)
    dist = canvas.scaleunit * canvas.mbwidthf
    seg = trimbond(canvas, Segment{Point2D}(canvas.coords, u, v), uvis, vvis)
    seg1 = translate(seg, pi / 2, dist)
    seg2 = translate(seg, -pi / 2, dist)
    drawline!(canvas, seg, ucolor, vcolor)
    drawline!(canvas, seg1, ucolor, vcolor)
    drawline!(canvas, seg2, ucolor, vcolor)
    return
end

const BOND_DRAWER = Dict(
    1 => Dict(
        :none => singlebond!,
        :up => wedged!,
        :down => dashedwedged!,
        :revup => reversed_wedged!,
        :revdown => reversed_dashedwedged!,
        :unspecified => wavesingle!
    ),
    2 => Dict(
        :none => doublebond!,
        :clockwise => clockwisedouble!,
        :anticlockwise => counterdouble!,
        :cis_trans => crossdouble!
    ),
    3 => Dict(
        :none => triplebond!
    )
)

setbond!(
    canvas::Canvas, order, notation, u, v, ucolor, vcolor, uvis, vvis
) = BOND_DRAWER[order][notation](canvas, u, v, ucolor, vcolor, uvis, vvis)


"""
    draw2d!(canvas::Canvas, mol::UndirectedGraph; kwargs...)

Draw molecular image to the canvas.
"""
function draw2d!(canvas::Canvas, mol::SimpleMolGraph; kwargs...)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    # get coords
    if has_coords(mol)  # SDFAtom or has coordgen! precache
        crds = coords2d(mol)
        default_bond_style = (has_cache(mol, :e_coordgen_bond_style) ? 
            get_cache(mol, :e_coordgen_bond_style) : sdf_bond_style(mol))
    else  # default SMILESAtom
        crds, default_bond_style = coordgen(mol)
    end
    isempty(crds) && return
    # Canvas settings
    initcanvas!(canvas, crds, boundary(mol, crds))
    # Properties
    atomsymbol_ = atom_symbol(mol)
    charge_ = charge(mol)
    implicith_ = implicit_hydrogens(mol)
    bondorder_ = bond_order(mol)
    atomcolor_ = atom_color(mol; kwargs...)
    isatomvisible_ = is_atom_visible(mol; kwargs...)
    bondstyle_ = bond_style(mol.graph, bondorder_, default_bond_style, crds, sssr(mol))

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
                mk = atommarkupleft(canvas, atomsymbol_[i], charge_[i], implicith_[i])
                drawtextleft!(canvas, pos, mk, atomcolor_[i])
                continue
            elseif maximum(cosnbrs) < 0
                # >[atom]
                mk = atommarkupright(canvas, atomsymbol_[i], charge_[i], implicith_[i])
                drawtextright!(canvas, pos, mk, atomcolor_[i])
                continue
            end
        end
        # -[atom]- or no hydrogens
        mk = atommarkupright(canvas, atomsymbol_[i], charge_[i], implicith_[i])
        drawtextcenter!(canvas, pos, mk, atomcolor_[i])
    end
    return
end


function drawatomindex!(canvas::Canvas, isatomvisible, color, bgcolor)
    for (i, v) in enumerate(isatomvisible)
        offset = v ? (0.0, canvas.fontsize/2.0) : (0.0, 0.0)
        pos = Point2D(canvas.coords, i) + offset
        drawtextannot!(canvas, pos, string(i), color, bgcolor)
    end
end


function sethighlight!(canvas::Canvas, edge_list::Vector{<:Edge}, color)
    for e in edge_list
        seg = Segment{Point2D}(canvas.coords, src(e), dst(e))
        drawlinehighlight!(canvas, seg, color)
    end
end

function sethighlight!(canvas::Canvas, node_list::Vector{<:Integer}, color)
    for i in node_list
        pos = Point2D(canvas.coords, i)
        drawtexthighlight!(canvas, pos, color)
    end
end
