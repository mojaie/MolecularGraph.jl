#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


function optpos(g::SimpleGraph, coords::Vector{Point2d})
    oppos = Vector{Point2d}(undef, length(coords))
    opdeg = Vector{Float64}(undef, length(coords))
    for i in vertices(g)
        oppos[i], opdeg[i] = optpos(coords[i], [coords[n] for n in neighbors(g, i)])
    end
    return oppos, opdeg
end


function optpos(center::Point2d, coords::Vector{Point2d})
    if length(coords) == 0
        return (Point2d(NaN, NaN), 0.0)
    elseif length(coords) == 1
        return (normalize(coords[1] - center) * -1, 0.0)
    elseif length(coords) == 2
        return (normalize(normalize(coords[1] - center) + normalize(coords[2] - center)) * -1, 0.0)
    else
        angs = [atan(reverse(coords[n] - center)...) for n in 1:length(coords)]
        sp = sortperm(angs, rev=true)
        diffs = Float64[]
        for k in 1:(length(coords) - 1)
            push!(diffs, angs[sp[k]] - angs[sp[k + 1]])
        end
        push!(diffs, angs[sp[end]] - angs[sp[1]] + 2 * pi)
        dmax, mi = findmax(diffs)
        a, b = mi == length(coords) ? (mi, 1) : (mi, mi + 1)
        mid = normalize(normalize(coords[a] - center) + normalize(coords[b] - center))
        return (dmax < pi ? mid : mid * -1, dmax)
    end
end


function label_direction(
        g::SimpleGraph, opos::Vector{Point2d}, odeg::Vector{Float64})
    arr = Vector{Symbol}(undef, length(odeg))
    for i in vertices(g)
        if isnan(opos[i][1])  # isolated node(ex. H2O, HCl)
            arr[i] = :right
        elseif (degree(g, i) > 2 && odeg[i] < pi) || (degree(g, i) == 2 && abs(opos[i][2]) > 0.87)  # -[atom]-
            arr[i] = :center
        elseif opos[i][1] < 0  # [atom]<
            arr[i] = :left
        else # >[atom]
            arr[i] = :right
        end
    end
    return arr
end


"""
    is_atom_visible(mol::SimpleMolGraph; show_carbon=:simple, kwargs...) -> Vector{Bool}

Return whether the atom is visible in the 2D drawing.
`show_carbon`: `simple` - does not show carbon labels, `terminal` - show only terminal carbonlabels,
`all` - show all carbon labels.
"""
is_atom_visible(
    mol::SimpleMolGraph; show_carbon=:simple, kwargs...
) = is_atom_visible(
    mol.graph, atom_symbol(mol), atom_charge(mol), multiplicity(mol),
    isotope(mol), bond_order(mol); kwargs...
)

function is_atom_visible(
        g::SimpleGraph, sym::Vector{Symbol}, chg::Vector{Int}, mul::Vector{Int},
        iso::Vector{Int}, bo::Vector{Int}; show_carbon=:simple, kwargs...)
    arr = trues(nv(g))
    show_carbon === :all && return arr
    ernk = edge_rank(g)
    for i in vertices(g)
        sym[i] === :C || continue
        chg[i] == 0 || continue
        mul[i] == 1 || continue
        iso[i] == 0 || continue
        degree(g, i) == 0 && continue
        degree(g, i) == 1 && show_carbon === :terminal && continue
        if degree(g, i) == 2
            nbrs = neighbors(g, i)
            u = edge_rank(ernk, i, nbrs[1])
            v = edge_rank(ernk, i, nbrs[2])
            if (bo[u] == 2 && bo[v] == 2)
                continue # allene-like
            end
        end
        arr[i] = false
    end
    return arr
end


function atom_markup_ordered(
        atommarkup::Vector{Vector{Vector{Tuple{Symbol,String}}}},
        labeldir::Vector{Symbol}, chg::Vector{Int})
    arr = Vector{Vector{Tuple{Symbol,String}}}(undef, length(chg))
    for i in 1:length(chg)
        if labeldir[i] === :left
            arr[i] = vcat(reverse(atommarkup[i])...)
        else
            arr[i] = vcat(atommarkup[i]...)
        end
        if chg[i] != 0
            push!(arr[i], (:sup, chargesign(chg[i])))
        end
    end
    return arr
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


function double_bond_style(
        g::SimpleGraph, bondorder_::Vector{Int}, coords::Vector{Point2d},
        sssr_::Vector{Vector{Int}})
    arr = Vector{Symbol}(undef, ne(g))
    ernk = edge_rank(g)
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
            se = edge_rank(ernk, src(e), snbr)
            bondorder_[se] == 2
        end
        ddbs = map(dnbrs) do dnbr
            de = edge_rank(ernk, dst(e), dnbr)
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
            e = edge_rank(ernk, rr[i], rr[i + 1])
            bondorder_[e] == 2 || continue
            arr[e] = rr[i] < rr[i + 1] ? :clockwise : :anticlockwise
        end
    end
    return arr
end


function bond_style(
        g::SimpleGraph, bondorder::Vector{Int}, defaultbondstyle::Vector{Symbol},
        coords_::Vector{Point2d}, sssr_::Vector{Vector{Int}})
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
    normalize_coords(
        g::SimpleGraph, coords::Vector{Point2d},
        label_length::Vector{Float64}, label_direction::Vector{Symbol},
        fontsizef::Float64, scalef::Float64, paddingXf::Float64, paddingYf::Float64)
    ) -> (coords, width, height)

Get boundaries and an appropriate bond length unit for the molecule drawing
canvas.
"""
function normalize_coords(
        g::SimpleGraph, coords::Vector{Point2d},
        label_length::Vector{Int}, label_direction::Vector{Symbol},
        atomvisible::BitVector,
        fontsizef::Float64, scalef::Float64, paddingX::Float64, paddingY::Float64)
    # scalef: bond length in pixel
    # paddingX, paddingY: drawing area padding in pixel
    # fontsizef: estimated font width in pixel
    (left, right) = extrema([p[1] for p in coords])
    (bottom, top) = extrema([p[2] for p in coords])
    width = right - left
    height = top - bottom
    dists = []
    # Size unit
    for e in edges(g)
        d = norm(coords[dst(e)] - coords[src(e)])
        if d > 0.0001  # Remove overlapped
            push!(dists, d)
        end
    end
    if isempty(dists)
        long = max(width, height)
        unit = long > 0.0001 ? long / sqrt(nv(g)) : 1
    else
        unit = Statistics.median(dists) # Median bond length
    end
    # extra offset for long text vertices
    sf = scalef / unit
    exleft = 0.0  # extra left offset
    exright = 0.0  # extra right offset
    for i in vertices(g)
        atomvisible[i] || continue
        label_length[i] > 1 || continue
        hf = label_direction[i] === :center ? 0.5 : 1
        xoff = label_length[i] * fontsizef * hf / sf
        if label_direction[i] !== :right && abs(left - coords[i][1]) < xoff
            exleft = max(exleft, xoff)
        end
        if label_direction[i] !== :left && abs(right - coords[i][1]) < xoff
            exright = max(exright, xoff)
        end
    end
    # transform
    conv = p -> (p - Point2d(left - exleft, top)) * Point2d(1, -1) * sf + Point2d(paddingX, paddingY)
    new_coords = conv.(coords)
    width = (width + exleft + exright) * sf + paddingX * 2
    height = height * sf + paddingY * 2
    return (new_coords, width, height)
end



"""
    draw2d!(canvas::Canvas, mol::UndirectedGraph; kwargs...)

Draw molecular image to the canvas.
"""
function draw2d!(canvas::Canvas, mol::SimpleMolGraph; kwargs...)
    # Calculate text direction, normalize coordinates and determine canvas size
    has_coords2d(mol) || coordgen!(mol)
    crds = coords2d(mol)
    oppos, opdeg = optpos(mol.graph, crds)
    labeldir = label_direction(mol.graph, oppos, opdeg)
    atommarkup = atom_markup_ordered(atom_markup(mol), labeldir, atom_charge(mol))
    charcnt = [length(join(e[2] for e in amk)) for amk in atommarkup]
    isatomvisible_ = is_atom_visible(mol; kwargs...)
    ncrds, width, height = normalize_coords(
        mol.graph, crds, charcnt, labeldir, isatomvisible_,
        canvas.fontsize * 0.5, canvas.scaleunit, canvas.paddingX, canvas.paddingY)
    initcanvas!(canvas, ncrds, width, height)
    # Draw bonds first, and then draw atoms
    bondorder_ = bond_order(mol)
    atomcolor_ = atom_color(mol; kwargs...)
    bondstyle_ = bond_style(mol.graph, bondorder_, draw2d_bond_style(mol), crds, sssr(mol))
    for (i, e) in enumerate(edges(mol))
        setbond!(
            canvas, bondorder_[i], bondstyle_[i], e.src, e.dst,
            atomcolor_[e.src], atomcolor_[e.dst],
            isatomvisible_[e.src], isatomvisible_[e.dst]
        )
    end
    for i in vertices(mol)
        isatomvisible_[i] || continue
        alabel = atom_label(atommarkup[i], canvas.fonttagmap)
        drawtext!(canvas, canvas.coords[i], alabel, atomcolor_[i], labeldir[i])
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
