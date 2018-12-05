#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


# TODO: deprecated

export
    default_annot!


struct IntrinsicAnnot <: Annotation end
struct ValenceAnnot <: Annotation end
struct AtomAnnot <: Annotation end
struct GroupAnnot <: Annotation end


function fragment_annot!(mol)
    required_annotation(mol, :GroupAnnot)

    # Carboxylic acid
    # 1: C, 2: =O, 3: -OH(or -[O-])
    # Ester
    # 1: C, 2: =O, 3: -O-
    mol.v[:Carboxyl] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    mol.v[:Ester] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for nbrs in mol.graph.adjacency[mol.v[:Carbonyl] .== 2]
        ns = collect(keys(nbrs))
        if mol.v[:Oxygen] == 4 in ns
            if mol.v[:Oxygen] == 1 in ns
                mol.v[:Carboxyl][n] =
            elseif mol.v[:Oxygen] == 2 in ns
                mol.v[:Ester][n] =
            end
        end
    end
    # Peroxide
    # 1: -O-, 2: -OH
    mol.v[:Peroxide] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for nbrs in mol.graph.adjacency[mol.v[:Oxygen] .== 1 or 2]
        ns = collect(keys(nbrs))
        if n1 mol.v[:Symbol] == :O
            if mol.v[:Oxygen] == 1 in ns
                mol.v[:Peroxide][n] = 2
            elseif mol.v[:Oxygen] == 2 in ns
                mol.v[:Peroxide][n] = 1
            end
        end
    end
    # Acetal
    # 1: C, 2: -OR, 3: -OH
    mol.v[:Peroxide] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for nbrs in mol.graph.adjacency[mol.v[:Oxygen] .== 1 or 2]
        ns = collect(keys(nbrs))
        if n1 mol.v[:Symbol] == :O
            if mol.v[:Oxygen] == 1 in ns
                mol.v[:Peroxide][n] = 2
            elseif mol.v[:Oxygen] == 2 in ns
                mol.v[:Peroxide][n] = 1
            end
        end
    end


end
