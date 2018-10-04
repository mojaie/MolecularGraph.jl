#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

function assignvalence(mol::MolecularGraph)
    for (u, v, bond) in bondsiter(mol)
        uobj = getatom(mol, u)
        vobj = getatom(mol, v)
        if bond.order == 2
            uobj.pi = 1
            vobj.pi = 1
            vobj.carbonylC = uobj.symbol == "O" && !uobj.charge
            uobj.carbonylC = vobj.symbol == "O" && !vobj.charge
        elseif bond.order == 3
            uobj.pi = vobj.pi = 2
        end
    end
    maxnbrs = Dict([
        ("C", 4), ("Si", 4), ("N", 3), ("P", 3), ("As", 3),
        ("O", 2), ("S", 2), ("Se", 2), ("F", 1), ("Cl", 1), ("Br", 1), ("I", 1)
    ])
    for (i, nbrs) in neighborsiter(mol)
        iobj = getatom(mol, i)
        if length(nbrs) == 2 && all(bond.order == 2 for bond in nbrs.values()
            iobj.pi = 2 # sp (allene, ketene)
        end
        if iobj.symbol in maxnbrs
            Hs = maxnbrs[iobj.symbol] - length(nbrs) - iobj.pi + iobj.charge
            if Hs > 0
                addhydrogen!(iobj, Hs)
            end
        end
    end
    add(mol.descriptors, "Valence")
end


function assign_rotatable(mol::MolecularGraph)
    required_descriptor(mol, "Valence")
    required_descriptor(mol, "Topology")
end
