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


function wedgeduv!(canvas, uv, ucolor, vcolor)
    drawwedge!(canvas, swap(uv), ucolor)
    return
end


function wedgedvu!(canvas, uv, ucolor, vcolor)
    drawwedge!(canvas, uv, ucolor)
    return
end


function dashedwedgeduv!(canvas, uv, ucolor, vcolor)
    drawdashedwedge!(canvas, swap(uv), ucolor)
    return
end


function dashedwedgedvu!(canvas, uv, ucolor, vcolor)
    drawdashedwedge!(canvas, uv, ucolor)
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


function ringdouble!(canvas::Canvas, uv, ucolor, vcolor, rad)
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


function triplebond!(canvas::Canvas, uv, ucolor, vcolor)
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
        upos = vec2d(coords[atompos(mol, bond.u), :])
        vpos = vec2d(coords[atompos(mol, bond.v), :])
        uv = vecpair(upos, vpos)
        if upos == vpos
            continue # avoid zero division
        end
        u = uatom.visible ? vecU(trimU(uv, canvas.trimoverlapf)) : upos
        v = vatom.visible ? vecV(trimV(uv, canvas.trimoverlapf)) : vpos
        drawer = BOND_DRAWER[bond.order][bond.notation]
        drawer(canvas, vecpair(u, v),
               Color(getcolor(uatom)...), Color(getcolor(vatom)...))
    end
    """ Draw atoms """
    for (i, atom) in enumerate(atomvector(mol))
        if !atom.visible
            continue
        end
        pos = vec2d(coords[i, :])
        # Determine text direction
        if atom.Hcount > 0
            cosnbrs = []
            hrzn = vec2d(posX(pos) + 1, posY(pos))
            for nbr in keys(neighbors(mol, atom.index))
                posnbr = vec2d(coords[atompos(mol, nbr), :])
                dist = norm(posnbr - pos)
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
