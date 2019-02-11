#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    draw2d!,
    drawatomindex!


"""
    draw2d!(canvas::Canvas, mol::VectorMol;
            setting=copy(DRAW_SETTING), recalculate=false)

Draw molecular image to the canvas.
"""
function draw2d!(canvas::Canvas, mol::VectorMol;
                 setting=copy(DRAW_SETTING), recalculate=false)
    elemental!(mol)
    topology!(mol)

    # 2D coordinate generation
    if recalculate
        mol.coords[:Cartesian2D] = coords2d(mol)
    elseif !haskey(mol.coords, :Cartesian2D)
        if nodetype(mol) === SmilesAtom
            mol.coords[:Cartesian2D] = coords2d(mol)
        else
            matrix = zeros(Float64, nodecount(mol), 2)
            for (i, node) in nodesiter(mol)
                matrix[i, :] = node.coords[1:2]
            end
            mol.coords[:Cartesian2D] = cartesian2d(matrix)
        end
    end

    # Canvas settings
    atomnotation2d!(mol, setting=setting)
    bondnotation2d!(mol, setting=setting)
    initcanvas!(canvas, mol)
    if !canvas.valid
        return
    end

    # Draw bonds
    for (i, bond) in edgesiter(mol)
        u = bond.u
        v = bond.v
        drawer = BOND_DRAWER[mol[:BondOrder][i]][mol[:BondNotation][i]]
        drawer(canvas, mol, u, v)
    end

    # Draw atoms
    for i in 1:atomcount(mol)
        if !mol[:Visible2D][i]
            continue
        end
        pos = _point(mol.coords[:Cartesian2D], i)
        # Determine text direction
        if mol[:H_Count][i] > 0
            cosnbrs = []
            hrzn = [pos[1] + 1.0, pos[2]]
            for nbr in neighborkeys(mol, i)
                posnbr = _point(mol.coords[:Cartesian2D], nbr)
                dist = norm(posnbr - pos)
                if dist > 0
                    dp = dot(hrzn - pos, posnbr - pos)
                    push!(cosnbrs, dp / dist)
                end
            end
            if isempty(cosnbrs) || minimum(cosnbrs) > 0
                # [atom]< or isolated node(ex. H2O, HCl)
                atomsymbolright!(canvas, mol, i)
                continue
            elseif maximum(cosnbrs) < 0
                # >[atom]
                atomsymbolleft!(canvas, mol, i)
                continue
            end
        end
        # -[atom]- or no hydrogens
        atomsymbolcenter!(canvas, mol, i)
    end
    return
end


function drawatomindex!(canvas::Canvas, mol::VectorMol;
                        color=Color(0, 0, 0), bgcolor=Color(240, 240, 255))
    for i in 1:atomcount(mol)
        offset = mol[:Visible2D][i] ? [0 canvas.fontsize/2] : [0 0]
        pos = _point(mol.coords[:Cartesian2D], i) + offset
        atom_annotation!(canvas, pos, i, color, bgcolor)
    end
    return
end
