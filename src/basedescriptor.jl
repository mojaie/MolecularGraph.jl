#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    assign_descriptors!,
    assign_valence!,
    assign_rotatable!,
    rotatable_count


function assign_descriptors!(mol::MolecularGraph)
    molgraph_topology!(mol)
    assign_valence!(mol)
    assign_rotatable!(mol)
    return
end


struct Valence <: Descriptor end
struct Rotatable <: Descriptor end


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
    mol.descriptor[:Valence] = Valence()
    return
end


function assign_rotatable!(mol::MolecularGraph)
    required_descriptor(mol, :Valence)
    required_descriptor(mol, :Topology)
    ringmap = mol.descriptor[:Topology].cyclemap
    for b in bondvector(mol)
        isec = length(intersect(ringmap[b.u], ringmap[b.v]))
        ulen = length(neighbors(mol, b.u))
        vlen = length(neighbors(mol, b.v))
        if b.order == 1 && isec == 0 && ulen > 1 && vlen > 1
            b.rotatable = true
        end
    end
    mol.descriptor[:Rotatable] = Rotatable()
end


function rotatable_count(mol::MolecularGraph)
    reduce(+, 1 for b in bondvector(mol) if b.rotatable; init=0)
end
