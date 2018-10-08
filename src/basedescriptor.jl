#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

function assign_valence(mol::MolecularGraph)
    for bond in mol.graph.bonds
        u = getatom(mol, bond.u)
        v = getatom(mol, bond.v)
        if bond.order == 2
            u.pi = 1
            v.pi = 1
            v.carbonylC = u.symbol == "O" && !u.charge
            u.carbonylC = v.symbol == "O" && !v.charge
        elseif bond.order == 3
            u.pi = v.pi = 2
        end
    end
    maxnbrs = Dict([
        ("C", 4), ("Si", 4), ("N", 3), ("P", 3), ("As", 3),
        ("O", 2), ("S", 2), ("Se", 2), ("F", 1), ("Cl", 1), ("Br", 1), ("I", 1)
    ])
    for (i, nbrs) in mol.graph.adjmap
        atom = getatom(mol, i)
        if length(nbrs) == 2 && all(b.order == 2 for b in values(nbrs))
            atom.pi = 2 # sp (allene, ketene)
        end
        if atom.symbol in maxnbrs
            Hs = maxnbrs[atom.symbol] - length(nbrs) - atom.pi + atom.charge
            if Hs > 0
                addhydrogen!(atom, Hs)
            end
        end
    end
    add(mol.descriptors, "Valence")
end
