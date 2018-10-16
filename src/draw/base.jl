#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

using Statistics

export
    Canvas,
    Color,
    readytodraw!,
    display_terminal_carbon!,
    equalize_terminal_double_bond!,
    double_bond_along_ring!,
    coordsvector,
    booundary,
    sizeunit


abstract type Canvas end


struct Color
    r::UInt8
    g::UInt8
    b::UInt8
end


function readytodraw!(mol::MolecularGraph)
    required_descriptor(mol, "Valence")
    required_descriptor(mol, "Topology")
    # display_terminal_carbon!(mol)
    equalize_terminal_double_bond!(mol)
    double_bond_along_ring!(mol)
end


function display_terminal_carbon!(mol::MolecularGraph)
    for atom in atomvector(mol)
        if length(neighbors(mol, atom.index)) == 1
            atom.visible = true
        end
    end
    return
end


function equalize_terminal_double_bond!(mol::MolecularGraph)
    for atom in atomvector(mol)
        nbrs = neighbors(mol, atom.index)
        if length(nbrs) == 1
            nbr = collect(values(nbrs))[1]
            if nbr.order == 2
                nbr.notation = 2
            end
        end
    end
    return
end


function double_bond_along_ring!(mol::MolecularGraph)
    for ring in sort(mol.rings, by=length, rev=true)
        vtcs = [point2d(getatom(mol, n).coords[1:2]) for n in ring.arr]
        cw = isclockwise(vtcs)
        if cw === nothing
            continue
        end
        directed = cw ? ring.arr : reverse(ring.arr)
        push!(directed, directed[1])
        for i in 1: length(directed) - 1
            u = directed[i]
            v = directed[i + 1]
            bond = getbond(mol, u, v)
            if bond.order == 2
                bond.notation = u > v ? 0 : 1
            end
        end
    end
    return
end


function coordsvector(mol::MolecularGraph)
    atoms = atomvector(mol)
    coords = zeros(Float32, length(atoms), 2)
    for (i, atom) in enumerate(atoms)
        coords[i, :] = collect(atom.coords[1:2])
    end
    coords
end


function boundary(coords::Matrix{Float32})
    (left, right) = extrema(coords[:, 1])
    (bottom, top) = extrema(coords[:, 2])
    width = right - left
    height = top - bottom
    (top, left, width, height)
end


function sizeunit(mol::MolecularGraph, coords::Matrix{Float32})
    dists = []
    for bond in bondvector(mol)
        upos = point2d(coords[atompos(mol, bond.u), :])
        vpos = point2d(coords[atompos(mol, bond.v), :])
        d = distance(upos, vpos)
        if d > 0  # Remove overlapped
            push!(dists, d)
        end
    end
    # Median bond length
    isempty(dists) ? nothing : median(dists)
end
