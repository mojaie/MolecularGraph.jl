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


function wedgeduv!(canvas, seg, ucolor, vcolor)
    drawwedge!(canvas, Segment(seg.v, seg.u), ucolor)
    return
end


function wedgedvu!(canvas, seg, ucolor, vcolor)
    drawwedge!(canvas, seg, ucolor)
    return
end


function dashedwedgeduv!(canvas, seg, ucolor, vcolor)
    drawdashedwedge!(canvas, Segment(seg.v, seg.u), ucolor)
    return
end


function dashedwedgedvu!(canvas, seg, ucolor, vcolor)
    drawdashedwedge!(canvas, seg, ucolor)
    return
end


function wavesingle!(canvas, seg, ucolor, vcolor)
    drawwave!(canvas, seg, ucolor)
    return
end


function doublebond!(canvas, seg, ucolor, vcolor)
    dist = canvas.scalef * canvas.mbwidthf / 2
    seg1 = translate(seg, pi / 2, dist)
    seg2 = translate(seg, -pi / 2, dist)
    drawline!(canvas, seg1, ucolor, vcolor)
    drawline!(canvas, seg2, ucolor, vcolor)
    return
end


function crossdouble!(canvas, seg, ucolor, vcolor)
    dist = canvas.scalef * canvas.mbwidthf / 2
    seg1 = translate(seg, pi / 2, dist)
    seg2 = translate(seg, -pi / 2, dist)
    drawline!(canvas, Segment(seg1.u, seg2.v), ucolor, vcolor)
    drawline!(canvas, Segment(seg2.u, seg1.v), ucolor, vcolor)
    return
end


function ringdouble!(canvas::Canvas, seg, ucolor, vcolor, rad)
    dist = canvas.scalef * canvas.mbwidthf
    segin = translate(seg, rad, dist)
    segtr = trim_uv(segin, canvas.triminnerf)
    drawline!(canvas, seg, ucolor, vcolor)
    drawline!(canvas, segtr, ucolor, vcolor)
    return
end

clockwisedouble!(canvas, seg, ucolor, vcolor) = ringdouble!(
    canvas, seg, ucolor, vcolor, -pi / 2)

counterdouble!(canvas, seg, ucolor, vcolor) = ringdouble!(
    canvas, seg, ucolor, vcolor, pi / 2)


function triplebond!(canvas::Canvas, seg, ucolor, vcolor)
    dist = canvas.scalef * canvas.mbwidthf
    seg1 = translate(seg, pi / 2, dist)
    seg2 = translate(seg, -pi / 2, dist)
    drawline!(canvas, seg, ucolor, vcolor)
    drawline!(canvas, seg1, ucolor, vcolor)
    drawline!(canvas, seg2, ucolor, vcolor)
    return
end


BOND_DRAWER = Dict(
    1 => Dict(
        0 => singlebond!,
        1 => wedgeduv!,
        2 => wedgedvu!,
        3 => dashedwedgeduv!,
        4 => dashedwedgedvu!,
        5 => wavesingle!
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
                drawtext!(canvas, pos, atom, :right)
                continue
            elseif maximum(cosnbrs) < 0
                # >[atom]
                drawtext!(canvas, pos, atom, :left)
                continue
            end
        end
        # -[atom]- or no hydrogens
        drawtext!(canvas, pos, atom, :center)
    end
end
