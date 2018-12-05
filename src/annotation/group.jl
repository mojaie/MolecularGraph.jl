#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export group_annot!


struct GroupAnnot <: Annotation end


function group_annot!(mol)
    required_annotation(mol, :AtomType)

    # Alchohol
    # 0: methanol, 1: primary, 2: sec, 3:tert, 4: vinyl
    mol.v[:Alcohol] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for nbrs in mol.graph.adjacency[mol.v[:Oxygen] .== 1]
        n = collect(keys(nbrs))[1]
        if mol.v[:Symbol][n] != :C
            continue
        elseif mol.v[:Pi][n] == 1
            mol.v[:Alcohol][i] = 4
        else
            mol.v[:Alcohol][n] = mol.v[:Degree][n] - 1
        end
    end

    # Thiol
    # 0: methanethiol, 1: primary, 2: sec, 3:tert, 4: vinyl
    mol.v[:Thiol] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for nbrs in mol.graph.adjacency[mol.v[:Sulfur] .== 1]
        n = collect(keys(nbrs))[1]
        if mol.v[:Symbol][n] != :C
            continue
        elseif mol.v[:Pi][n] == 1
            mol.v[:Thiol][i] = 4
        else
            mol.v[:Thiol][n] = mol.v[:Degree][n] - 1
        end
    end

    # Carbonyl
    # 0:formaldehyde, 1: aldehyde, 2: ketone
    mol.v[:Carbonyl] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for nbrs in mol.graph.adjacency[mol.v[:Oxygen] .== 4]
        n = collect(keys(nbrs))[1]
        if mol.v[:Symbol][n] != :C
            continue
        else
            mol.v[:Carbonyl][n] = mol.v[:Degree][n] - 1
        end
    end

    # Imine
    # 0:methaneimine, 1:aldimine, 2:ketimine
    mol.v[:Imine] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for nbrs in mol.graph.adjacency[mol.v[:Nitrogen] .== 5]
        n = [i for (i, b) in nbrs if mol.v[:BondOrder][b] == 2][1]
        if mol.v[:Symbol][n] != :C
            continue
        else
            mol.v[:Imine][n] = mol.v[:Degree][n] - 1
        end
    end

    # Sulfoxide
    # 1:Sulfoxide 2:Sulfone
    mol.v[:Sulfoxide] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for nbrs in mol.graph.adjacency[mol.v[:Oxygen] .== 4]
        n = collect(keys(nbrs))[1]
        if mol.v[:Symbol][n] != :S
            continue
        else
            mol.v[:Sulfoxide][n] = mol.v[:Sulfoxide][n] === nothing ? 1 : 2
        end
    end

    mol.annotation[:AtomGroup] = GroupAnnot()
end
