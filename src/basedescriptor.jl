#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    assign_descriptors!,
    assign_valence!


function assign_descriptors!(mol::MolecularGraph)
    assign_valence!(mol)
    topology!(mol)
    minifyring!(mol)
    return
end


function assign_valence!(mol::MolecularGraph)
    for bond in bondvector(mol)
        u = getatom(mol, bond.u)
        v = getatom(mol, bond.v)
        if bond.order == 2
            u.pi = 1
            v.pi = 1
            v.carbonylC = u.symbol == "O" && u.charge == 0
            u.carbonylC = v.symbol == "O" && v.charge == 0
        elseif bond.order == 3
            u.pi = v.pi = 2
        end
    end
    maxnbrs = Dict([
        ("C", 4), ("Si", 4), ("N", 3), ("P", 3), ("As", 3),
        ("O", 2), ("S", 2), ("Se", 2), ("F", 1), ("Cl", 1), ("Br", 1), ("I", 1)
    ])
    for (i, nbrs) in enumerate(adjvector(mol))
        atom = mol.graph.nodes[i]
        if length(nbrs) == 2 && all(b.order == 2 for b in values(nbrs))
            atom.pi = 2 # sp (allene, ketene)
        end
        if atom.symbol in keys(maxnbrs)
            Hs = maxnbrs[atom.symbol] - length(nbrs) - atom.pi + atom.charge
            if Hs > 0
                sethydrogen!(atom, convert(UInt8, Hs))
            end
        end
    end
    push!(mol.descriptors, "Valence")
    return
end
