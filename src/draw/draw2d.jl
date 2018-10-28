#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    draw!


function singlebond!(canvas, uv, ucolor, vcolor)
    drawline!(canvas, uv, ucolor, vcolor)
    return
end


function wedged!(canvas, uv, ucolor, vcolor)
    drawwedge!(canvas, swap(uv), ucolor)
    return
end


function dashedwedged!(canvas, uv, ucolor, vcolor)
    drawdashedwedge!(canvas, swap(uv), ucolor)
    return
end


function wavesingle!(canvas, uv, ucolor, vcolor)
    drawwave!(canvas, uv, ucolor)
    return
end


function doublebond!(canvas, uv, ucolor, vcolor)
    dist = canvas.scalef * canvas.mbwidthf / 2
    uv1 = translate(uv, pi / 2, dist)
    uv2 = translate(uv, -pi / 2, dist)
    drawline!(canvas, uv1, ucolor, vcolor)
    drawline!(canvas, uv2, ucolor, vcolor)
    return
end


function crossdouble!(canvas, uv, ucolor, vcolor)
    dist = canvas.scalef * canvas.mbwidthf / 2
    uv1 = translate(uv, pi / 2, dist)
    uv2 = translate(uv, -pi / 2, dist)
    drawline!(canvas, vecpair(vecU(uv1), vecV(uv2)), ucolor, vcolor)
    drawline!(canvas, vecpair(vecU(uv2), vecV(uv1)), ucolor, vcolor)
    return
end


function ringdouble!(canvas, uv, ucolor, vcolor, rad)
    dist = canvas.scalef * canvas.mbwidthf
    uvin = translate(uv, rad, dist)
    uvtr = trimUV(uvin, canvas.triminnerf)
    drawline!(canvas, uv, ucolor, vcolor)
    drawline!(canvas, uvtr, ucolor, vcolor)
    return
end

clockwisedouble!(canvas, uv, ucolor, vcolor) = ringdouble!(
    canvas, uv, ucolor, vcolor, -pi / 2)

counterdouble!(canvas, uv, ucolor, vcolor) = ringdouble!(
    canvas, uv, ucolor, vcolor, pi / 2)


function triplebond!(canvas, uv, ucolor, vcolor)
    dist = canvas.scalef * canvas.mbwidthf
    uv1 = translate(uv, pi / 2, dist)
    uv2 = translate(uv, -pi / 2, dist)
    drawline!(canvas, uv, ucolor, vcolor)
    drawline!(canvas, uv1, ucolor, vcolor)
    drawline!(canvas, uv2, ucolor, vcolor)
    return
end


BOND_DRAWER = Dict(
    1 => Dict(
        0 => singlebond!,
        1 => wedged!,
        6 => dashedwedged!,
        4 => wavesingle!
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


function draw!(canvas::Canvas, mol::Molecule)
    if atomcount(mol) == 0
        return
    end
    initialize!(canvas, mol)
    coords = canvas.coords

    """ Draw bonds """
    for (i, bond) in enumerate(mol.graph.edges)
        upos = vec2d(coords[bond.u, :])
        vpos = vec2d(coords[bond.v, :])
        uv = vecpair(upos, vpos)
        if upos == vpos
            continue # Overlapped: avoid zero division
        end
        u = mol.v[:AtomVisible][bond.u] ? vecU(trimU(uv, canvas.trimoverlapf)) : upos
        v = mol.v[:AtomVisible][bond.v] ? vecV(trimV(uv, canvas.trimoverlapf)) : vpos
        drawer = BOND_DRAWER[mol.v[:BondOrder][i]][mol.v[:BondNotation][i]]
        drawer(canvas, vecpair(u, v),
               mol.v[:AtomColor][bond.u], mol.v[:AtomColor][bond.v])
    end
    """ Draw atoms """
    for i in 1:atomcount(mol)
        if !mol.v[:AtomVisible][i]
            continue
        end
        pos = vec2d(coords[i, :])
        # Determine text direction
        if mol.v[:H_Count][i] > 0
            cosnbrs = []
            hrzn = vec2d(posX(pos) + 1, posY(pos))
            for nbr in keys(neighbors(mol, i))
                posnbr = vec2d(coords[nbr, :])
                dist = norm(posnbr - pos)
                if dist > 0
                    dp = dot(hrzn - pos, posnbr - pos)
                    push!(cosnbrs, dp / dist)
                end
            end
            if isempty(cosnbrs) || minimum(cosnbrs) > 0
                # [atom]< or isolated node(ex. H2O, HCl)
                drawrighttext!(
                    canvas, pos, mol.v[:Symbol][i], mol.v[:Charge][i],
                    mol.v[:H_Count][i], mol.v[:AtomColor][i]
                )
                continue
            elseif maximum(cosnbrs) < 0
                # >[atom]
                drawlefttext!(
                    canvas, pos, mol.v[:Symbol][i], mol.v[:Charge][i],
                    mol.v[:H_Count][i], mol.v[:AtomColor][i]
                )
                continue
            end
        end
        # -[atom]- or no hydrogens
        drawcentertext!(
            canvas, pos, mol.v[:Symbol][i], mol.v[:Charge][i],
            mol.v[:H_Count][i], mol.v[:AtomColor][i]
        )
    end
end
