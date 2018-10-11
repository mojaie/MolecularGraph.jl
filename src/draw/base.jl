#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


export
    Canvas,
    draw!


abstract type Canvas end

CANVAS = Dict(:svg => SvgCanvas)

function draw!(mol::MolecularGraph, format::Symbol)
    required_descriptor(mol, "Valence")
    required_descriptor(mol, "Topology")
    # display_terminal_carbon!(mol)
    equalize_terminal_double_bond!(mol)
    double_bond_along_ring!(mol)
    canvas = CANVAS[format]()
    draw!(canvas, mol)
    show(canvas)
end


# TODO: non-destructive draw (clone mol)


function display_terminal_carbon!(mol::MolecularGraph)
    adj = adjmap(mol)
    for atom in atomvector(mol)
        if length(adj[atom.index]) == 1
            atom.visible = true
        end
    end
    return
end


function equalize_terminal_double_bond!(mol::MolecularGraph)
    adj = adjmap(mol)
    for atom in atomvector(mol)
        if length(adj[atom.index]) == 1
            nbr = values(adj[atom.index])[1]
            if nbr.order == 2
                nbr.notation = 2
            end
        end
    end
    return
end


function double_bond_along_ring!(mol::MolecularGraph)
    for ring in sort(mol.rings, by=length, rev=true)
        vtcs = [Point2D(getatom(mol, n).coords)[1:2] for n in ring.arr]
        cw = isclockwise(vtcs)
        if cw === nothing
            continue
        directed = cw ? ring.arr : reverse(ring.arr)
        push!(directed, directed[1])
        for i in 1: length(directed)
            bond = getbond(mol, directed[i], directed[i + 1])
            if bond.order == 2
                bond.notation = u > v ? 1 : 2
            end
        end
    end
    return
end
