#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

BOND_TRIM = 0.3

BOND_DRAWER = {
    1: {
        0: singlebond,
        1: wedgedsingle,
        2: dashedwedgedsingle,
        3: wavesingle,
    },
    2: {
        0: clockwisedouble,
        1: counterdouble,
        2: doublebond,
        3: crossdouble
    },
    3: {
        0: triplebond
    }
}


function draw!(canvas::Canvas, mol::MolecularGraph)
    scaleandcenter!(canvas, mol)
    mlb = mol.size2d[3] # median bond length
    if length(atomvector(mol)) == 0
        return
    end

    """ Draw bonds """
    for bond in bondvector(mol)
        if !bond.visible
            continue
        end
        uatom = getatom(mol, bond.u)
        vatom = getatom(mol, bond.v)
        upos = uatom.coords
        vpos = vatom.coords
        if upos == vpos
            continue # avoid zero division
        end
        if uatom.visible
            upos = paralleltrim(upos, vpos, BOND_TRIM, 2)[1]
        end
        if vatom.visible
            vpos = paralleltrim(upos, vpos, BOND_TRIM, 1)[2]
        end
        color = uatom.color
        vcolor = vatom.color
        drawer = BOND_DRAWER[bond.order][bond.notation]
        drawer(canvas, p1, p2, color, vcolor=vcolor, mlb)
    end

    """ Draw atoms """
    for atom in atomvedtor(mol)
        if !atom.visible
            continue
        end
        pos = atom.coords
        color = atom.color
        # Determine text direction
        if atom.Hcount > 0
            cosnbrs = []
            hrzn = (pos[1] + 1, pos[2])
            for nbr in keys(neighbors(mol, atom.index))
                posnbr = getatom(mol, nbr).coords
                dist = distance(pos, posnbr)
                if dist > 0
                    dp = dot(hrzn, posnbr, pos)
                    push!(cosnbrs, dp / dist)
                end
            end
            if isempty(cosnbrs) || minimum(cosnbrs) > 0
                # [atom]< or isolated node(ex. H2O, HCl)
                text = atom.htmlformula(:right)
                drawtext(canvas, pos, text, color, :right)
                continue
            elseif maximum(cosnbrs) < 0
                # >[atom]
                text = atom.htmlformula(:left)
                drawtext(canvas, pos, text, color, :left)
                continue
            end
        end
        # -[atom]- or no hydrogens
        text = atom.htmlformula(:left)
        drawtext(canvas, pos, text, color, :center)
    end
end


function scaleandcenter!(canvas::Canvas, mol::MolecularGraph)
    atomcount = length(atomvector(mol))
    if atomcount < 2
        canvas.size2d = (0, 0, 1)
        return
    end
    xs = []
    ys = []
    for atom in atomvector(mol)
        push!(xs, atom.coords[1])
        push!(ys, atom.coords[2])



function singlebond(canvas::Canvas, p1, p2, color, vcolor, mlb)
    drawline(canvas, u, v, color, vcolor)
end
