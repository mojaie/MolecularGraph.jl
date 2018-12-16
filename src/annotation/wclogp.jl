#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    WCLogP,
    wclogp!,
    wclogptype!,
    wclogpcontrib!


const WCLOGP_TABLE = YAML.load(open(
    joinpath(dirname(@__FILE__), "..", "..", "assets", "const", "wclogp.yaml")
))


struct WCLogP <: Annotation end


function wclogp!(mol)
    required_annotation(mol, :FuncGroup)
    wclogptype!(mol)
    wclogpcontrib!(mol)
    mol.annotation[:WCLogP] = WCLogP()
    return
end


function wclogptype!(mol)
    mol.v[:WCLogP] = zeros(atomcount(mol))

    # Aliphatic carbon (C1-4,8-12)
    for i in findall(mol.v[:Symbol] .== :C .& mol.v[:Pi] .== 0)
        nbrs = neighborkeys(mol, i)
        if all(mol.v[:Aromatic][nbrs] .== false)
            if all(mol.v[:Symbol][nbrs] .== :C)
                if mol.v[:Degree][i] <= 2
                    mol.v[:WCLogP][i] = :C1
                else
                    mol.v[:WCLogP][i] = :C2
                end
            else
                if mol.v[:Degree][i] <= 2
                    mol.v[:WCLogP][i] = :C3
                else
                    mol.v[:WCLogP][i] = :C4
                end
            end
        else
            if mol.v[:Degree][i] == 1
                nbr = nbrs[1]
                if mol.v[:Symbol][nbr] == :C
                    mol.v[:WCLogP][i] = :C8
                else
                    mol.v[:WCLogP][i] = :C9
                end
            elseif mol.v[:Degree][i] == 2
                mol.v[:WCLogP][i] = :C10
            elseif mol.v[:Degree][i] == 3
                mol.v[:WCLogP][i] = :C11
            elseif mol.v[:Degree][i] == 4
                mol.v[:WCLogP][i] = :C12
            end
        end
    end
    # Aliphatic multiple bond carbon (C5-7,26)
    for i in findall(
            mol.v[:Symbol] .== :C .* mol.v[:Pi] .== 1 .* mol.v[:Aromatic] .== false)
        if :Alkyne in keys(mol.v) && mol.v[:Alkyne] == 1
            mol.v[:WCLogP][i] = :C7
            continue
        elseif :Nitrile in keys(mol.v) && mol.v[:Nitrile] == 1
            mol.v[:WCLogP][i] = :C7
            continue
        elseif :Alkene in keys(mol.v) && mol.v[:Alkene] == 1
            nbrs = neighborkeys(mol, i)
            if all(mol.v[:Aromatic][nbrs] .== false)
                mol.v[:WCLogP][i] = :C6
            else
                mol.v[:WCLogP][i] = :C26
            end
        end
    end
    # Aromatic carbon (C14-25)
    for i in findall(mol.v[:Symbol] .== :C .* mol.v[:Aromatic] .== true)

    end

    # Exotic carbon (C13,27,S)

    # Nitrogen (N1-14,S)

end


function wclogpcontrib!(mol)
    mol.v[:WCLogPContrib] = zeros(atomcount(mol))
    for i in atomkeys
        hcnt = mol.v[:H_Count][i]
        if hcnt == 0
            continue
        end
        contrib = WCLOGP_TABLE[string(mol.v[:WCLogP][i])]
        if mol.v[:Symbol][i] == :C
            # H1
        elseif :Carboxyl in keys(mol.v) && mol.v[:Carboxyl][i] == 3
            # H4
        elseif :Alcohol in keys(mol.v) && mol.v[:Alcohol][i] == 2
            # H2
        elseif :Amine in keys(mol.v) && mol.v[:Amine][i] == 2
            # H3
        else
            # HS
        mol.v[:WCLogPContrib] = contrib
    end
end

"""
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
"""
