#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Aromatic,
    aromatic!


struct Aromatic <: Annotation end


function aromatic!(mol)
    required_annotation(mol, :Topology)
    required_annotation(mol, :Elemental)
    # Precalculate carbonyl
    mol.v[:Aromatic] = falses(atomcount(mol))
    mol.v[:AromaticBond] = falses(bondcount(mol))
    for ring in mol.annotation[:Topology].rings
        sub = nodesubgraph(mol.graph, Set(ring))
        if satisfyHuckel(mol, ring)
            mol.v[:Aromatic][ring] .= true
            for e in edgekeys(sub)
                mol.v[:AromaticBond][e] = true
            end
        elseif mol isa GVectorMol{SmilesAtom,SmilesBond}
            # SMILES aromatic atom
            for (i, n) in nodesiter(sub)
                if n.isaromatic
                    mol.v[:Aromatic][i] = true
                end
            end
            for (i, e) in edgesiter(sub)
                if e.isaromatic
                    mol.v[:AromaticBond][i] = true
                end
            end
        end
    end
    mol.annotation[:Aromatic] = Aromatic()
    return
end


function satisfyHuckel(mol::VectorMol, ring)
    cnt = 0
    carbonyl = fgroupquery(mol, "[#6]=[OD1]")
    if !isempty(carbonyl)
        carbonyl = union(carbonyl...)
    end
    for r in ring
        if r in carbonyl
            continue
        elseif mol.v[:Pi][r] == 1
            cnt += 1
        elseif mol.v[:LonePair][r] === nothing
            return false
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
