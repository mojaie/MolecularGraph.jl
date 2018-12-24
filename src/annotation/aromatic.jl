#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Aromatic,
    aromatic!


struct Aromatic <: Annotation end


function aromatic!(mol::VectorMol)
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
    # Carbonyl
    carbonylC = Int[]
    carbonylO = findall((mol.v[:Symbol] .== :O) .* (mol.v[:Degree] .== 1))
    for o in carbonylO
        c = collect(neighborkeys(mol.graph, o))[1]
        if mol.v[:Symbol][c] == :C
            push!(carbonylC, c)
        end
    end
    for r in ring
        if r in carbonylC
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
