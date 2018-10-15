#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


BOND_DRAWER = {
    1: {
        1: singlebond,
        2: wedgedsingle,
        3: dashedwedgedsingle,
        4: wavesingle,
    },
    2: {
        1: clockwisedouble,
        2: counterdouble,
        3: doublebond,
        4: crossdouble
    },
    3: {
        1: triplebond
    }
}


function draw!(canvas::Canvas, mol::MolecularGraph)
    if length(atomvector(mol)) == 0
        finalize!(canvas)
        return
    end
    # Scale and center
    (coords, width, height, unit) = scaleandcenter(mol)
    canvas.width = width
    canvas.height = height
    canvas.unit = unit

    """ Draw bonds """
    for bond in bondvector(mol)
        if !bond.visible
            continue
        end
        uatom = getatom(mol, bond.u)
        vatom = getatom(mol, bond.v)
        upos = coords[atompos(mol, bond.u), :]
        vpos = coords[atompos(mol, bond.v), :]
        if upos == vpos
            continue # avoid zero division
        end
        u = uatom.visible
            ? trim_u(Segment(upos, vpos), canvas.trimoverlapf)[1] : upos
        v = vatom.visible
            ? trim_v(Segment(upos, vpos), canvas.trimoverlapf)[2] : vpos
        drawer = BOND_DRAWER[bond.order][bond.notation]
        drawer(canvas, u, v, uatom.color, vcolor=vatom.color)
    end

    """ Draw atoms """
    for (i, atom) in enumerate(atomvector(mol))
        if !atom.visible
            continue
        end
        pos = coords[i, :]
        color = atom.color
        # Determine text direction
        if atom.Hcount > 0
            cosnbrs = []
            hrzn = (pos[1] + 1, pos[2])
            for nbr in keys(neighbors(mol, atom.index))
                posnbr = coords[atompos(mol, nbr), :]
                dist = distance(pos, posnbr)
                if dist > 0
                    dp = dot(hrzn, posnbr, pos)
                    push!(cosnbrs, dp / dist)
                end
            end
            if isempty(cosnbrs) || minimum(cosnbrs) > 0
                # [atom]< or isolated node(ex. H2O, HCl)
                text = atom.htmlformula(:right)
                drawtext! (canvas, pos, text, color, :right)
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
    finalize!(canvas)
end


function scaleandcenter(mol::MolecularGraph)
    atomcount = length(atomvector(mol))
    coords = zeros(Float32, atomcount, 2)
    if atomcount == 1
        return (coords, 0, 0, 1)
    end
    for (i, atom) in enumerate(atomvector(mol))
        coords[i, :] = collect(atom.coords[1:2])
    end
    (xmin, xmax) = extrema(coords[:, 1])
    (ymin, ymax) = extrema(coords[:, 2])
    width = xmax - xmin
    height = ymax - ymin
    xoffset = width / 2 + xmin
    yoffset = height / 2 + ymin
    dists = []
    for bond in bondvector(mol)
        upos = Point2D(coords[atompos(mol, bond.u), :])
        vpos = Point2D(coords[atompos(mol, bond.v), :])
        d = distance(upos, vpos)
        if d > 0  # Remove overlapped
            push!(dists, d)
        end
    end
    # Median bond length (as a drawing size unit)
    medianbondlength = isempty(dists)
        ? sqrt(max(width, height) / atomcount) : median(dists)
    # Centering
    centered = broadcast(-, coords, [xoffset yoffset])

    return (centered, width, height, medianbondlength)
end


function singlebond!(canvas::Canvas, u, v, color, vcolor)
    drawline!(canvas, u, v, color, vcolor)
    return
end


function wedgedsingle!(canvas::Canvas, u, v, color, vcolor)
    drawwedge!(canvas, u, v, color)
    return
end


function dashedwedgedsingle!(canvas::Canvas, u, v, color, vcolor)
    drawdashedwedge!(canvas, u, v, color)
    return
end


function wavesingle!(canvas::Canvas, u, v, color, vcolor)
    drawwave!(canvas, u, v, color)
    return
end


function doublebond!(canvas::Canvas, u, v, color, vcolor)
    dist = canvas.mbwidthf / 2 * canvas.unit
    (u1, v1) = trim_uv_move(Segment(u, v), true, dist, 1)
    (u2, v2) = trim_uv_move(Segment(u, v), false, dist, 1)
    drawline!(canvas, u1, v1, color, vcolor)
    drawline!(canvas, u2, v2, color, vcolor)
    return
end


function crossdouble!(canvas::Canvas, u, v, color, vcolor)
    dist = canvas.mbwidthf / 2 * canvas.unit
    (u1, v1) = trim_uv_move(Segment(u, v), true, dist, 1)
    (u2, v2) = trim_uv_move(Segment(u, v), false, dist, 1)
    drawline!(canvas, u1, v2, color, vcolor)
    drawline!(canvas, u2, v1, color, vcolor)
    return
end


function ringdouble!(canvas::Canvas, u, v, color, vcolor, direction)
    dist = canvas.mbwidthf * canvas.unit
    (u1, v1) = trim_uv_move(Segment(u, v), direction, dist, canvas.triminnerf)
    drawline!(canvas, u, v, color, vcolor)
    drawline!(canvas, u1, v1, color, vcolor)
    return
end

clockwisedouble!(canvas, u, v, color, vcolor) = ringdouble(
    canvas, u, v, color, vcolor, true)

counterdouble!(canvas, u, v, color, vcolor) = ringdouble(
    canvas, u, v, color, vcolor, false)


function triplebond!(canvas::Canvas, u, v, color, vcolor)
    dist = canvas.mbwidthf * canvas.unit
    (u1, v1) = trim_uv_move(Segment(u, v), true, dist, 1)
    (u2, v2) = trim_uv_move(Segment(u, v), false, dist, 1)
    drawline!(canvas, u, v, color, vcolor)
    drawline!(canvas, u1, v1, color, vcolor)
    drawline!(canvas, u2, v2, color, vcolor)
    return
end
