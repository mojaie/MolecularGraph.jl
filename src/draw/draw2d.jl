#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

"""
    is_atom_visible(mol::SimpleMolGraph; show_carbon=:simple, kwargs...) -> Vector{Bool}

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
    dispatch_update!(mol)
    return is_atom_visible(mol.graph, atom_symbol(mol), atom_charge(mol), multiplicity(mol),
        [atom_mass(props(mol, i)) for i in vertices(mol)], bond_order(mol); kwargs...)
end


"""
    atom_label(
        contents::Vector{Tuple{Symbol,String}},
        mapping::Dict{Symbol,Tuple{String,String}}
    ) -> String

Generate atom label from ordered and annotated texts and the markup tag mapping specific to
the drawing canvas.
"""
function atom_label(
        contents::Vector{Tuple{Symbol,String}},
        mapping::Dict{Symbol,Tuple{String,String}})
    texts = String[]
    for (tag, cont) in contents
        if tag === :default
            push!(texts, cont)
        else
            st, ed = mapping[tag]
            push!(texts, st, cont, ed)
        end
    end
    return join(texts)
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
        cw = isclockwise(coords[ring])
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
    boundary(mol::SimpleMolGraph, coords::AbstractArray{Float64}
        ) -> (top, left, width, height, unit)

Get boundaries and an appropriate bond length unit for the molecule drawing
canvas.
"""
function boundary(mol::SimpleMolGraph, coords::Vector{Point2d})
    (left, right) = extrema([p[1] for p in coords])
    (bottom, top) = extrema([p[2] for p in coords])
    width = right - left
    height = top - bottom
    dists = []
    # Size unit
    for e in edges(mol)
        d = norm(coords[dst(e)] - coords[src(e)])
        if d > 0.0001  # Remove overlapped
            push!(dists, d)
        end
    end
    if isempty(dists)
        long = max(width, height)
        unit = long > 0.0001 ? long / sqrt(nv(mol)) : 1
    else
        unit = Statistics.median(dists) # Median bond length
    end
    return (top, left, width, height, unit)
end


function trimbond_(uvis, vvis)
    uvis && vvis && return trim_uv
    uvis && return trim_u
    vvis && return trim_v
    return (seg, t) -> seg
end

function singlebond!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis)
    seg = trimbond_(uvis, vvis)(Line(canvas.coords[[u, v]]...), canvas.trimoverlap)
    drawline!(canvas, seg, ucolor, vcolor)
    return
end

function wedged!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis)
    seg = trimbond_(uvis, vvis)(Line(canvas.coords[[u, v]]...), canvas.trimoverlap)
    drawwedge!(canvas, seg, ucolor, vcolor)
    return
end

function reversed_wedged!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis)
    seg = trimbond_(vvis, uvis)(Line(canvas.coords[[v, u]]...), canvas.trimoverlap)
    drawwedge!(canvas, seg, vcolor, ucolor)
    return
end

function dashedwedged!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis)
    seg = trimbond_(uvis, vvis)(Line(canvas.coords[[u, v]]...), canvas.trimoverlap)
    drawdashedwedge!(canvas, seg, ucolor, vcolor)
    return
end

function reversed_dashedwedged!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis)
    seg = trimbond_(vvis, uvis)(Line(canvas.coords[[v, u]]...), canvas.trimoverlap)
    drawdashedwedge!(canvas, seg, vcolor, ucolor)
    return
end

function wavesingle!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis)
    seg = trimbond_(uvis, vvis)(Line(canvas.coords[[u, v]]...), canvas.trimoverlap)
    drawwave!(canvas, seg, ucolor, vcolor)
    return
end

function doublebond!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis)
    seg = trimbond_(uvis, vvis)(Line(canvas.coords[[u, v]]...), canvas.trimoverlap)
    seg1 = translate(seg, pi / 2, canvas.mbwidth / 2)
    seg2 = translate(seg, -pi / 2, canvas.mbwidth / 2)
    drawline!(canvas, seg1, ucolor, vcolor)
    drawline!(canvas, seg2, ucolor, vcolor)
    return
end

function crossdouble!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis)
    seg = trimbond_(uvis, vvis)(Line(canvas.coords[[u, v]]...), canvas.trimoverlap)
    cw = translate(seg, pi / 2, canvas.mbwidth / 2)
    ccw = translate(seg, -pi / 2, canvas.mbwidth / 2)
    drawline!(canvas, Line(cw[1], ccw[2]), ucolor, vcolor)
    drawline!(canvas, Line(ccw[1], cw[2]), ucolor, vcolor)
    return
end

function ringdouble!(canvas::Canvas, u, v, ucolor, vcolor, uvis, vvis, direction)
    seg = trimbond_(uvis, vvis)(Line(canvas.coords[[u, v]]...), canvas.trimoverlap)
    segin = trimbond_(!uvis, !vvis)(
        translate(seg, direction, canvas.mbwidth), canvas.triminner)
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
    seg = trimbond_(uvis, vvis)(Line(canvas.coords[[u, v]]...), canvas.trimoverlap)
    seg1 = translate(seg, pi / 2, canvas.mbwidth)
    seg2 = translate(seg, -pi / 2, canvas.mbwidth)
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
    dispatch_update!(mol)
    # 2D coordinates required
    if !has_coords2d(mol)
        coordgen!(mol)
    end
    crds = coords2d(mol)
    # Canvas settings
    initcanvas!(canvas, crds, boundary(mol, crds))
    # Properties
    implicith_ = implicit_hydrogens(mol)
    bondorder_ = bond_order(mol)
    atomcolor_ = atom_color(mol; kwargs...)
    isatomvisible_ = is_atom_visible(mol; kwargs...)
    bondstyle_ = bond_style(mol.graph, bondorder_, draw2d_bond_style(mol), crds, sssr(mol))
    atommarkup_ = atom_markup(mol)

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
        pos = canvas.coords[i]
        # Determine text direction
        if implicith_[i] > 0
            cosnbrs = []
            hrzn = pos + Point2d(1, 0)
            for nbr in neighbors(mol, i)
                posnbr = canvas.coords[nbr]
                dist = norm(posnbr - pos)
                if dist > 0
                    dp = dot(hrzn - pos, posnbr - pos)
                    push!(cosnbrs, dp / dist)
                end
            end
            if isempty(cosnbrs) || minimum(cosnbrs) > 0
                # [atom]< or isolated node(ex. H2O, HCl)
                alabel = atom_label(vcat(reverse(atommarkup_[i])...), canvas.fonttagmap)
                drawtext!(canvas, pos, alabel, atomcolor_[i], :left)
                continue
            elseif maximum(cosnbrs) < 0
                # >[atom]
                alabel = atom_label(vcat(atommarkup_[i]...), canvas.fonttagmap)
                drawtext!(canvas, pos, alabel, atomcolor_[i], :right)
                continue
            end
        end
        # -[atom]- or no hydrogens
        alabel = atom_label(vcat(atommarkup_[i]...), canvas.fonttagmap)
        drawtext!(canvas, pos, alabel, atomcolor_[i], :center)
    end
    return
end


function drawatomindex!(canvas::Canvas, isatomvisible, color, bgcolor)
    for (i, v) in enumerate(isatomvisible)
        offset = v ? Point2d(0, canvas.fontsize/2.0) : Point2d(0, 0)
        pos = canvas.coords[i] + offset
        drawtextannot!(canvas, pos, string(i), color, bgcolor)
    end
end


function sethighlight!(canvas::Canvas, edge_list::Vector{<:Edge}, color)
    for e in edge_list
        seg = Line(canvas.coords[[src(e), dst(e)]]...)
        drawlinehighlight!(canvas, seg, color)
    end
end

function sethighlight!(canvas::Canvas, node_list::Vector{<:Integer}, color)
    for i in node_list
        pos = canvas.coords[i]
        drawtexthighlight!(canvas, pos, color)
    end
end
