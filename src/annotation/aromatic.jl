#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export aromatic!


function aromatic!(mol)
    required_annotation(mol, :Topology)
    required_annotation(mol, :AtomGroup)
    mol.v[:Aromatic] = falses(atomcount(mol))
    for ring in mol.annotation[:Topology].rings
        if satisfyHuckel(mol, ring)
            mol.v[:Aromatic][ring] .= true
        end
    end
    return
end


function satisfyHuckel(mol::VectorMol, ring)
    cnt = 0
    for r in ring
        if mol.v[:Carbonyl][r] == 2
            continue
        elseif mol.v[:Pi][r] == 1
            cnt += 1
        elseif mol.v[:LonePair][r] > 0
            cnt += 2
        elseif mol.v[:LonePair][r] < 0
            continue
        else
            return false
        end
    end
    cnt % 4 == 2
end
