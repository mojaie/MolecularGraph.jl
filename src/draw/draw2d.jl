#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    draw!


function singlebond!(canvas::Canvas, segment, color, vcolor)
    drawline!(canvas, segment, color, vcolor)
    return
end


function wedgedsingle!(canvas::Canvas, seg, color, vcolor)
    drawwedge!(canvas, seg, color)
    return
end


function dashedwedgedsingle!(canvas::Canvas, seg, color, vcolor)
    drawdashedwedge!(canvas, seg, color)
    return
end


function wavesingle!(canvas::Canvas, seg, color, vcolor)
    drawwave!(canvas, seg, color)
    return
end


function doublebond!(canvas::Canvas, seg, color, vcolor)
    dist = canvas.mbwidthf / 2 * length(seg)
    seg1 = trim_uv_move(seg, true, dist, 1)
    seg2 = trim_uv_move(seg, false, dist, 1)
    drawline!(canvas, seg1, color, vcolor)
    drawline!(canvas, seg2, color, vcolor)
    return
end


function crossdouble!(canvas::Canvas, seg, color, vcolor)
    dist = canvas.mbwidthf / 2 * length(seg)
    seg1 = trim_uv_move(seg, true, dist, 1)
    seg2 = trim_uv_move(seg, false, dist, 1)
    drawline!(canvas, Segment(seg1.u, seg2.v), color, vcolor)
    drawline!(canvas, Segment(seg2.u, seg1.v), color, vcolor)
    return
end


function ringdouble!(canvas::Canvas, seg, color, vcolor, direction)
    dist = canvas.mbwidthf * length(seg)
    segin = trim_uv_move(seg, direction, dist, canvas.triminnerf)
    drawline!(canvas, seg, color, vcolor)
    drawline!(canvas, segin, color, vcolor)
    return
end

clockwisedouble!(canvas, seg, color, vcolor) = ringdouble!(
    canvas, seg, color, vcolor, true)

counterdouble!(canvas, seg, color, vcolor) = ringdouble!(
    canvas, seg, color, vcolor, false)


function triplebond!(canvas::Canvas, seg, color, vcolor)
    dist = canvas.mbwidthf * length(seg)
    seg1 = trim_uv_move(seg, true, dist, 1)
    seg2 = trim_uv_move(seg, false, dist, 1)
    drawline!(canvas, seg, color, vcolor)
    drawline!(canvas, seg1, color, vcolor)
    drawline!(canvas, seg2, color, vcolor)
    return
end


BOND_DRAWER = Dict(
    1 => Dict(
        1 => singlebond!,
        2 => wedgedsingle!,
        3 => dashedwedgedsingle!,
        4 => wavesingle!
    ),
    2 => Dict(
        1 => clockwisedouble!,
        2 => counterdouble!,
        3 => doublebond!,
        4 => crossdouble!
    ),
    3 => Dict(
        1 => triplebond!
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
        u = uatom.visible ? trim_u(Segment(upos, vpos), canvas.trimoverlapf)[1] : upos
        v = vatom.visible ? trim_v(Segment(upos, vpos), canvas.trimoverlapf)[2] : vpos
        drawer = BOND_DRAWER[bond.order][bond.notation]
        drawer(canvas, Segment(u, v), getcolor(uatom), getcolor(vatom))
    end

    """ Draw atoms """
    for (i, atom) in enumerate(atomvector(mol))
        if !atom.visible
            continue
        end
        pos = point2d(coords[i, :])
        color = getcolor(atom)
        # Determine text direction
        if atom.Hcount > 0
            cosnbrs = []
            hrzn = (pos[1] + 1, pos[2])
            for nbr in keys(neighbors(mol, atom.index))
                posnbr = point2d(coords[atompos(mol, nbr), :])
                dist = distance(pos, posnbr)
                if dist > 0
                    dp = dot(hrzn, posnbr, pos)
                    push!(cosnbrs, dp / dist)
                end
            end
            if isempty(cosnbrs) || minimum(cosnbrs) > 0
                # [atom]< or isolated node(ex. H2O, HCl)
                text = atom.htmlformula(:right)
                drawtext!(canvas, pos, text, color, :right)
                continue
            elseif maximum(cosnbrs) < 0
                # >[atom]
                text = atom.htmlformula(:left)
                drawtext!(canvas, pos, text, color, :left)
                continue
            end
        end
        # -[atom]- or no hydrogens
        text = atom.htmlformula(:left)
        drawtext!(canvas, pos, text, color, :center)
    end
end
