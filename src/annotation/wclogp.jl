#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    WCLogP,
    wclogpcalc!,
    wclogptype!,
    wclogphydrogentype,
    wclogpcontrib!


const WCLOGP_TABLE = YAML.load(open(
    joinpath(dirname(@__FILE__), "..", "..", "assets", "const", "wclogp.yaml")
))["logP"]

struct WCLogP <: Annotation end


"""
    wclogpcalc!(mol::VectorMol)

Assign Wildman-Crippen LogP atom types and contributions.
"""
function wclogpcalc!(mol::VectorMol)
    haskey(mol, :WCLogPContrib) && return
    aromatic!(mol)
    wclogptype!(mol)
    wclogpcontrib!(mol)
    return
end


function wclogptype!(mol::VectorMol)
    mol[:WCLogP] = fill(:undef, atomcount(mol))
    ALIPH_HETERO = Set([:H, :C, :N, :O, :P, :S, :F, :Cl, :Br, :I])
    AROM_HETERO = Dict(
        :H  => :C18, :C => :C21, :N => :C22, :O => :C23, :S => :C24,
        :F => :C14, :Cl => :C15, :Br => :C16, :I => :C17
    )
    # TODO:
    P_BLOCK = (:Al, :B, :Si, :Ga, :Ge, :As, :Se, :Sn, :Te, :Pb,
               :Ne, :Ar, :Kr, :Xe, :Rn)
    D_BLOCK = (:Fe, :Co, :Cu, :Zn, :Tc, :Cd, :Pt, :Au, :Hg, :Gd)

    # Carbons
    for i in findall(mol[:Symbol] .== :C)
        if mol[:Pi][i] == 0
            # Aliphatic (C1-4,8-12,27)
            nbrs = collect(neighborkeys(mol, i))
            if !isempty(setdiff(mol[:Symbol][nbrs], ALIPH_HETERO))
                mol[:WCLogP][i] = :C27 # Adjacent to inorganic atoms
            elseif all(mol[:Aromatic][nbrs] .== false)
                if all(mol[:Symbol][nbrs] .== :C)
                    # Aliphatic carbon (C1,2)
                    if mol[:Degree][i] <= 2
                        mol[:WCLogP][i] = :C1
                    else
                        mol[:WCLogP][i] = :C2
                    end
                else
                    # Adjacent to heteroatoms (C3,4)
                    if mol[:Degree][i] <= 2
                        mol[:WCLogP][i] = :C3
                    else
                        mol[:WCLogP][i] = :C4
                    end
                end
            else
                # Adjacent to aromatic atoms (C8-12)
                if mol[:Degree][i] == 1
                    nbr = nbrs[1]
                    if mol[:Symbol][nbr] == :C
                        mol[:WCLogP][i] = :C8
                    else
                        mol[:WCLogP][i] = :C9
                    end
                elseif mol[:Degree][i] == 2
                    mol[:WCLogP][i] = :C10
                elseif mol[:Degree][i] == 3
                    mol[:WCLogP][i] = :C11
                elseif mol[:Degree][i] == 4
                    mol[:WCLogP][i] = :C12
                end
            end
        elseif mol[:Aromatic][i]
            # Aromatic (C13-25)
            if degree(mol, i) == 2
                mol[:WCLogP][i] = :C18 # Aromatic non-substituted
                continue
            end
            subst = nothing
            substbond = nothing
            aromcnt = 0
            for (nbr, b) in neighbors(mol, i)
                if mol[:AromaticBond][b]
                    aromcnt += 1
                else
                    subst = nbr
                    substbond = b
                end
            end
            if aromcnt == 3
                mol[:WCLogP][i] = :C19 # Bridgehead
            elseif !haskey(AROM_HETERO, mol[:Symbol][subst])
                mol[:WCLogP][i] = :C13 # Inorganic substituent
            elseif mol[:Aromatic][subst]
                mol[:WCLogP][i] = :C20 # Aromatic substituent
            elseif mol[:BondOrder][substbond] == 2
                mol[:WCLogP][i] = :C25 # Double bond substituent
            else
                # Typical substituent (C14-17,21-24)
                mol[:WCLogP][i] = AROM_HETERO[mol[:Symbol][subst]]
            end
        else
            # Aliphatic multiple bond (C5-7,26)
            nbrs = collect(neighborkeys(mol, i))
            bonds = collect(neighboredgekeys(mol, i))
            if any(mol[:BondOrder][bonds] .== 3)
                mol[:WCLogP][i] = :C7 # Alkyne, Nitrile
            elseif any((mol[:Pi][nbrs] .> 0) .* (mol[:Symbol][nbrs] .!= :C))
                mol[:WCLogP][i] = :C5 # Double bond to non-C
            elseif any(mol[:Aromatic][nbrs] .== true)
                mol[:WCLogP][i] = :C26 # Adjacent to aromatic
            else
                mol[:WCLogP][i] = :C6 # Double bond to C (including allene)
            end
        end
        if mol[:WCLogP][i] == :undef
            mol[:WCLogP][i] = :CS # Not found
        end
    end

    # Nitrogens
    for i in findall(mol[:Symbol] .== :N)
        if mol[:Charge][i] > 0
            nbrs = collect(neighborkeys(mol, i))
            if mol[:Aromatic][i]
                mol[:WCLogP][i] = :N12 # Charged aromatic nitrogen
            elseif mol[:H_Count][i] == 0
                if mol[:Degree][i] == 2 && all(mol[:Symbol][nbrs] .== :N)
                    mol[:WCLogP][i] = :N14 # Azide is exceptionally N14
                else
                    mol[:WCLogP][i] = :N13 # Quart-ammonium
                end
            elseif mol[:Pi][i] == 0
                mol[:WCLogP][i] = :N10 # Protonated amine
            else
                mol[:WCLogP][i] = :N14 # Other charged
            end
        elseif mol[:Charge][i] < 0
            mol[:WCLogP][i] = :N14 # Neg charged
        elseif mol[:Aromatic][i]
            mol[:WCLogP][i] = :N11 # Neutral aromatic
        elseif mol[:Pi][i] == 2
            mol[:WCLogP][i] = :N9 # Nitrile
        elseif mol[:Pi][i] == 1
            # Imine (N5,6)
            if mol[:Degree][i] == 1
                mol[:WCLogP][i] = :N5
            elseif mol[:Degree][i] == 2
                mol[:WCLogP][i] = :N6
            end
        else
            nbrs = collect(neighborkeys(mol, i))
            if all(mol[:Aromatic][nbrs] .== false)
                # Aliphatic amine (N1,2,7)
                if mol[:Degree][i] == 1
                    mol[:WCLogP][i] = :N1
                elseif mol[:Degree][i] == 2
                    mol[:WCLogP][i] = :N2
                elseif mol[:Degree][i] == 3
                    mol[:WCLogP][i] = :N7
                end
            else
                # Aromatic amine (N3,4,8)
                if mol[:Degree][i] == 1
                    mol[:WCLogP][i] = :N3
                elseif mol[:Degree][i] == 2
                    mol[:WCLogP][i] = :N4
                elseif mol[:Degree][i] == 3
                    mol[:WCLogP][i] = :N8
                end
            end
        end
        if mol[:WCLogP][i] == :undef
            # Not found
            mol[:WCLogP][i] = :NS
        end
    end

    # Oxygens
    for i in findall(mol[:Symbol] .== :O)
        if mol[:Aromatic][i]
            # Aromatic oxygen (O1)
            mol[:WCLogP][i] = :O1
            continue
        end
        nbrs = collect(neighborkeys(mol, i))
        bonds = collect(neighboredgekeys(mol, i))
        if mol[:Degree][i] == 1
            # Hydroxyl (O2,5-12)
            nbr = nbrs[1]
            bond = bonds[1]
            if mol[:WCLogP][nbr] == :C5 && mol[:Charge][i] < 0
                # Acid (O12)
                mol[:WCLogP][i] = :O12
            elseif mol[:BondOrder][bond] == 2 || mol[:Charge][i] < 0
                # Carbonyl or oxide (O5-11)
                if mol[:WCLogP][nbr] == :C25
                    mol[:WCLogP][i] = :O8 # Aromatic carbonyl
                elseif mol[:WCLogP][nbr] == :C5
                    # Aliphatic carbonyl
                    cnbrs = neighborkeys(mol, nbr)
                    pop!(cnbrs, i)
                    cnbrs = collect(cnbrs)
                    if all(mol[:Symbol][cnbrs] .!= :C)
                        mol[:WCLogP][i] = :O11 # Carbonyl heteroatom
                    elseif any(mol[:Aromatic][cnbrs] .== true)
                        mol[:WCLogP][i] = :O10 # Carbonyl aromatic
                    else
                        mol[:WCLogP][i] = :O9 # Carbonyl aliphatic
                    end
                elseif mol[:Symbol][nbr] in (:O, :N)
                    mol[:WCLogP][i] = :O5 # O2? or N-oxide
                elseif mol[:Symbol][nbr] == :S
                    mol[:WCLogP][i] = :O6 # S-oxide
                else
                    mol[:WCLogP][i] = :O7 # Other oxide
                end
            else
                # Alcohol (O2)
                mol[:WCLogP][i] = :O2
            end
        elseif mol[:Degree][i] == 2
            # Ether (O3,4)
            if all(mol[:Aromatic][nbrs] .== false)
                mol[:WCLogP][i] = :O3 # Aliphatic ether
            else
                mol[:WCLogP][i] = :O4 # Aromatic ether
            end
        end
    end

    # Others
    for i in findall(mol[:WCLogP] .== :undef)
        if mol[:Symbol][i] in (:F, :Cl, :Br, :I)
            if mol[:Charge][i] == 0
                mol[:WCLogP][i] = mol[:Symbol][i]
            else
                mol[:WCLogP][i] = :Hal
            end
        elseif mol[:Symbol][i] == :P
            mol[:WCLogP][i] = :P
        elseif mol[:Symbol][i] == :S
            if mol[:Aromatic][i]
                mol[:WCLogP][i] = :S3
            elseif mol[:Charge][i] == 0
                mol[:WCLogP][i] = :S1
            else
                mol[:WCLogP][i] = :S2
            end
        elseif mol[:Symbol][i] == :H
            nbr = pop!(neighborkeys(mol, i))
            mol[:WCLogP][i] = wclogphydrogentype(mol, nbr)
        elseif mol[:Symbol][i] in P_BLOCK
            mol[:WCLogP][i] = :Me1
        elseif mol[:Symbol][i] in D_BLOCK
            mol[:WCLogP][i] = :Me2
        end
    end
end


function wclogphydrogentype(mol::VectorMol, i)
    if mol[:Symbol][i] == :C
        return :H1 # Hydrocarbon
    elseif mol[:Symbol][i] == :N
        return :H3 # Amine
    elseif mol[:Symbol][i] == :O
        nbr = pop!(neighborkeys(mol, i))
        if mol[:Symbol][nbr] == :N
            return :H3 # Hydroxyamine
        elseif mol[:Symbol][nbr] in (:O, :S)
            return :H4 # Peroxide, sulfoxide
        elseif mol[:Symbol][nbr] == :C
            if mol[:Pi][nbr] == 1 && !mol[:Aromatic][nbr]
                return :H4 # Acid
            else
                return :H2 # Alcohol
            end
        else
            return :H2 # Other hydroxyl
        end
    else
        return :H2 # Other hydrate
    end
    return :HS # HS may be impossible
end


function wclogpcontrib!(mol::VectorMol)
    mol[:WCLogPContrib] = zeros(atomcount(mol))
    for i in nodekeys(mol)
        hcnt = mol[:H_Count][i]
        contrib = get(WCLOGP_TABLE, string(mol[:WCLogP][i]), 0)
        if mol[:Symbol][i] == :H
            mol[:WCLogPContrib][i] = 0  # avoid double count
        elseif hcnt == 0
            mol[:WCLogPContrib][i] = contrib
        else
            htype = wclogphydrogentype(mol, i)
            hcontrib = WCLOGP_TABLE[string(htype)] * hcnt
            mol[:WCLogPContrib][i] = contrib + hcontrib
        end
    end
    return
end
