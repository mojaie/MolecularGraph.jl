#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    draw!


function singlebond!(canvas, seg, ucolor, vcolor)
    drawline!(canvas, seg, ucolor, vcolor)
    return
end


function wedgedsingle!(canvas, seg, ucolor, vcolor)
    drawwedge!(canvas, seg, ucolor)
    return
end


function dashedwedgedsingle!(canvas, seg, ucolor, vcolor)
    drawdashedwedge!(canvas, seg, ucolor)
    return
end


function wavesingle!(canvas, seg, ucolor, vcolor)
    drawwave!(canvas, seg, ucolor)
    return
end


function doublebond!(canvas, seg, ucolor, vcolor)
    dist = canvas.mbwidthf / 2 * length(seg)
    seg1 = trim_uv_move(seg, true, dist, 1)
    seg2 = trim_uv_move(seg, false, dist, 1)
    drawline!(canvas, seg1, ucolor, vcolor)
    drawline!(canvas, seg2, ucolor, vcolor)
    return
end


function crossdouble!(canvas, seg, ucolor, vcolor)
    dist = canvas.mbwidthf / 2 * length(seg)
    seg1 = trim_uv_move(seg, true, dist, 1)
    seg2 = trim_uv_move(seg, false, dist, 1)
    drawline!(canvas, Segment(seg1.u, seg2.v), ucolor, vcolor)
    drawline!(canvas, Segment(seg2.u, seg1.v), ucolor, vcolor)
    return
end


function ringdouble!(canvas::Canvas, seg, ucolor, vcolor, direction)
    dist = canvas.mbwidthf * length(seg)
    segin = trim_uv_move(seg, direction, dist, canvas.triminnerf)
    drawline!(canvas, seg, ucolor, vcolor)
    drawline!(canvas, segin, ucolor, vcolor)
    return
end

clockwisedouble!(canvas, seg, ucolor, vcolor) = ringdouble!(
    canvas, seg, ucolor, vcolor, true)

counterdouble!(canvas, seg, ucolor, vcolor) = ringdouble!(
    canvas, seg, ucolor, vcolor, false)


function triplebond!(canvas::Canvas, seg, ucolor, vcolor)
    dist = canvas.mbwidthf * length(seg)
    seg1 = trim_uv_move(seg, true, dist, 1)
    seg2 = trim_uv_move(seg, false, dist, 1)
    drawline!(canvas, seg, ucolor, vcolor)
    drawline!(canvas, seg1, ucolor, vcolor)
    drawline!(canvas, seg2, ucolor, vcolor)
    return
end


BOND_DRAWER = Dict(
    1 => Dict(
        0 => singlebond!,
        1 => wedgedsingle!,
        2 => dashedwedgedsingle!,
        3 => wavesingle!
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


function draw!(canvas::Canvas, mol::MolecularGraph)
    if length(atomvector(mol)) == 0
        return
    end
    initialize!(canvas, mol)
    coords = canvas.coords

    """ Draw bonds """
    for bond in bondvector(mol)
        if !bond.visible
            continue
        end
        uatom = getatom(mol, bond.u)
        vatom = getatom(mol, bond.v)
        upos = point2d(coords[atompos(mol, bond.u), :])
        vpos = point2d(coords[atompos(mol, bond.v), :])
        if upos == vpos
            continue # avoid zero division
        end
        u = uatom.visible ? trim_u(Segment(upos, vpos), canvas.trimoverlapf).u : upos
        v = vatom.visible ? trim_v(Segment(upos, vpos), canvas.trimoverlapf).v : vpos
        drawer = BOND_DRAWER[bond.order][bond.notation]
        drawer(canvas, Segment(u, v),
               Color(getcolor(uatom)...), Color(getcolor(vatom)...))
    end

    """ Draw atoms """
    for (i, atom) in enumerate(atomvector(mol))
        if !atom.visible
            continue
        end
        pos = point2d(coords[i, :])
        color = Color(getcolor(atom)...)
        # Determine text direction
        if atom.Hcount > 0
            cosnbrs = []
            hrzn = Point2D(pos.x + 1, pos.y)
            for nbr in keys(neighbors(mol, atom.index))
                posnbr = point2d(coords[atompos(mol, nbr), :])
                dist = distance(pos, posnbr)
                if dist > 0
                    dp = dot(hrzn - pos, posnbr - pos)
                    push!(cosnbrs, dp / dist)
                end
            end
            if isempty(cosnbrs) || minimum(cosnbrs) > 0
                # [atom]< or isolated node(ex. H2O, HCl)
                text = htmlformula(atom, :right)
                drawtext!(canvas, pos, text, color, :right)
                continue
            elseif maximum(cosnbrs) < 0
                # >[atom]
                text = htmlformula(atom, :left)
                drawtext!(canvas, pos, text, color, :left)
                continue
            end
        end
        # -[atom]- or no hydrogens
        text = htmlformula(atom, :left)
        drawtext!(canvas, pos, text, color, :center)
    end
end
