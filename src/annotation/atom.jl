#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export atom_annot!


struct AtomAnnot <: Annotation end


function atom_annot!(mol)
    required_annotation(mol, :Default)

    # Oxygen type
    # 0:atom, 1:hydroxyl, 2:ether, 3:oxonium,
    # 4:oxo, 5:oxocarbenium
    mol.v[:Oxygen] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for i in 1:atomcount(mol)
        if mol.v[:Symbol][i] != :O
            continue
        end
        if mol.v[:Pi][i] == 1
            mol.v[:Oxygen][i] = mol.v[:Degree][i] + 3
        else
            mol.v[:Oxygen][i] = mol.v[:Degree][i]
        end
    end

    # Nitrogen type
    # 0:atom, 1:primary amine, 2:sec, 3:tert, 4:quart(ammonium),
    # 5:primary imine, 6:sec, 7:iminium, 8:nitrile
    mol.v[:Nitrogen] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for i in 1:atomcount(mol)
        if mol.v[:Symbol][i] != :N
            continue
        end
        if mol.v[:Pi][i] == 2
            mol.v[:Nitrogen][i] = 8
        elseif mol.v[:Pi][i] == 1
            mol.v[:Nitrogen][i] = mol.v[:Degree][i] + 4
        else
            mol.v[:Nitrogen][i] = mol.v[:Degree][i]
        end
    end

    # Sulfer type
    # 0:atom, 1:thiol, 2:sulfide, 3:sulfonium,
    # 4:thio, 5:thiocarbenium, 6:higher valent
    mol.v[:Sulfur] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for i in 1:atomcount(mol)
        if mol.v[:Symbol][i] != :N
            continue
        end
        if mol.v[:Valence][i] > 3
            mol.v[:Sulfur][i] = 6
        elseif mol.v[:Pi][i] == 1
            mol.v[:Sulfur][i] = mol.v[:Degree][i] + 4
        else
            mol.v[:Sulfur][i] = mol.v[:Degree][i]
        end
    end

    mol.annotation[:AtomType] = AtomAnnot()
end
