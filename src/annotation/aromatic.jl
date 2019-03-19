#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export aromatic!


"""
    aromatic!(mol::VectorMol)

Assign aromaticity parameters to the molecule parameter vector.

This sets two parameter vector `:Aromatic` and `:AromaticBond` to the molecule.
`:Aromatic` indicates whether the atoms are aromatic or not. Similarly,
`:AromaticBond` indicates whether the bonds are aromatic or not.

Note that aromaticity described here means simplified binary descriptor
(aromatic or not) based on classical Huckel's rule. This is intended for use in
some kind of pharmaceutical researches. Non-classical aromaticity such as
Moebius aromaticity is not considered.
"""
function aromatic!(mol::VectorMol)
    haskey(mol, :Aromatic) && return
    topology!(mol)
    elemental!(mol)
    # Precalculate carbonyl
    mol[:Aromatic] = falses(atomcount(mol))
    mol[:AromaticBond] = falses(bondcount(mol))
    for ring in mol.annotation[:Topology].rings
        sub = nodesubgraph(mol.graph, Set(ring))
        if satisfyHuckel(mol, ring)
            mol[:Aromatic][ring] .= true
            for e in edgekeys(sub)
                mol[:AromaticBond][e] = true
            end
        elseif nodetype(mol) === SmilesAtom
            # SMILES aromatic atom
            for (i, n) in nodesiter(sub)
                if n.isaromatic
                    mol[:Aromatic][i] = true
                end
            end
            for (i, e) in edgesiter(sub)
                if e.isaromatic
                    mol[:AromaticBond][i] = true
                end
            end
        end
    end
    return
end


function satisfyHuckel(mol::VectorMol, ring)
    cnt = 0
    # Carbonyl
    # TODO: polarized carbonyl
    carbonylC = Int[]
    carbonylO = findall(
        (mol[:Symbol] .== :O) .* (mol[:Degree] .== 1) .* (mol[:Pi] .== 1))
    for o in carbonylO
        c = collect(neighborset(mol.graph, o))[1]
        if mol[:Symbol][c] == :C
            push!(carbonylC, c)
        end
    end
    for r in ring
        if r in carbonylC
            continue
        elseif mol[:Pi][r] == 1
            cnt += 1
        elseif mol[:LonePair][r] === nothing
            return false
        elseif mol[:LonePair][r] > 0
            cnt += 2
        elseif mol[:LonePair][r] < 0
            continue
        else
            return false
        end
    end
    return cnt % 4 == 2
end
