#
# This file is a part of graphmol.jl
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


function wclogpcalc!(mol)
    required_annotation(mol, :Aromatic)
    wclogptype!(mol)
    wclogpcontrib!(mol)
    return
end


function wclogptype!(mol)
    mol.v[:WCLogP] = fill(:undef, atomcount(mol))
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
    for i in findall(mol.v[:Symbol] .== :C)
        if mol.v[:Pi][i] == 0
            # Aliphatic (C1-4,8-12,27)
            nbrs = collect(neighborkeys(mol, i))
            if !isempty(setdiff(mol.v[:Symbol][nbrs], ALIPH_HETERO))
                mol.v[:WCLogP][i] = :C27 # Adjacent to inorganic atoms
            elseif all(mol.v[:Aromatic][nbrs] .== false)
                if all(mol.v[:Symbol][nbrs] .== :C)
                    # Aliphatic carbon (C1,2)
                    if mol.v[:Degree][i] <= 2
                        mol.v[:WCLogP][i] = :C1
                    else
                        mol.v[:WCLogP][i] = :C2
                    end
                else
                    # Adjacent to heteroatoms (C3,4)
                    if mol.v[:Degree][i] <= 2
                        mol.v[:WCLogP][i] = :C3
                    else
                        mol.v[:WCLogP][i] = :C4
                    end
                end
            else
                # Adjacent to aromatic atoms (C8-12)
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
        elseif mol.v[:Aromatic][i]
            # Aromatic (C13-25)
            if degree(mol, i) == 2
                mol.v[:WCLogP][i] = :C18 # Aromatic non-substituted
                continue
            end
            subst = nothing
            substbond = nothing
            aromcnt = 0
            for (nbr, b) in neighbors(mol, i)
                if mol.v[:AromaticBond][b]
                    aromcnt += 1
                else
                    subst = nbr
                    substbond = b
                end
            end
            if aromcnt == 3
                mol.v[:WCLogP][i] = :C19 # Bridgehead
            elseif !haskey(AROM_HETERO, mol.v[:Symbol][subst])
                mol.v[:WCLogP][i] = :C13 # Inorganic substituent
            elseif mol.v[:Aromatic][subst]
                mol.v[:WCLogP][i] = :C20 # Aromatic substituent
            elseif mol.v[:BondOrder][substbond] == 2
                mol.v[:WCLogP][i] = :C25 # Double bond substituent
            else
                # Typical substituent (C14-17,21-24)
                mol.v[:WCLogP][i] = AROM_HETERO[mol.v[:Symbol][subst]]
            end
        else
            # Aliphatic multiple bond (C5-7,26)
            nbrs = collect(neighborkeys(mol, i))
            bonds = collect(neighboredgekeys(mol, i))
            if any(mol.v[:BondOrder][bonds] .== 3)
                mol.v[:WCLogP][i] = :C7 # Alkyne, Nitrile
            elseif any((mol.v[:Pi][nbrs] .> 0) .* (mol.v[:Symbol][nbrs] .!= :C))
                mol.v[:WCLogP][i] = :C5 # Double bond to non-C
            elseif any(mol.v[:Aromatic][nbrs] .== true)
                mol.v[:WCLogP][i] = :C26 # Adjacent to aromatic
            else
                mol.v[:WCLogP][i] = :C6 # Double bond to C (including allene)
            end
        end
        if mol.v[:WCLogP][i] == :undef
            mol.v[:WCLogP][i] = :CS # Not found
        end
    end

    # Nitrogens
    for i in findall(mol.v[:Symbol] .== :N)
        if mol.v[:Charge][i] > 0
            nbrs = collect(neighborkeys(mol, i))
            if mol.v[:Aromatic][i]
                mol.v[:WCLogP][i] = :N12 # Charged aromatic nitrogen
            elseif mol.v[:H_Count][i] == 0
                if mol.v[:Degree][i] == 2 && all(mol.v[:Symbol][nbrs] .== :N)
                    mol.v[:WCLogP][i] = :N14 # Azide is exceptionally N14
                else
                    mol.v[:WCLogP][i] = :N13 # Quart-ammonium
                end
            elseif mol.v[:Pi][i] == 0
                mol.v[:WCLogP][i] = :N10 # Protonated amine
            else
                mol.v[:WCLogP][i] = :N14 # Other charged
            end
        elseif mol.v[:Charge][i] < 0
            mol.v[:WCLogP][i] = :N14 # Neg charged
        elseif mol.v[:Aromatic][i]
            mol.v[:WCLogP][i] = :N11 # Neutral aromatic
        elseif mol.v[:Pi][i] == 2
            mol.v[:WCLogP][i] = :N9 # Nitrile
        elseif mol.v[:Pi][i] == 1
            # Imine (N5,6)
            if mol.v[:Degree][i] == 1
                mol.v[:WCLogP][i] = :N5
            elseif mol.v[:Degree][i] == 2
                mol.v[:WCLogP][i] = :N6
            end
        else
            nbrs = collect(neighborkeys(mol, i))
            if all(mol.v[:Aromatic][nbrs] .== false)
                # Aliphatic amine (N1,2,7)
                if mol.v[:Degree][i] == 1
                    mol.v[:WCLogP][i] = :N1
                elseif mol.v[:Degree][i] == 2
                    mol.v[:WCLogP][i] = :N2
                elseif mol.v[:Degree][i] == 3
                    mol.v[:WCLogP][i] = :N7
                end
            else
                # Aromatic amine (N3,4,8)
                if mol.v[:Degree][i] == 1
                    mol.v[:WCLogP][i] = :N3
                elseif mol.v[:Degree][i] == 2
                    mol.v[:WCLogP][i] = :N4
                elseif mol.v[:Degree][i] == 3
                    mol.v[:WCLogP][i] = :N8
                end
            end
        end
        if mol.v[:WCLogP][i] == :undef
            # Not found
            mol.v[:WCLogP][i] = :NS
        end
    end

    # Oxygens
    for i in findall(mol.v[:Symbol] .== :O)
        if mol.v[:Aromatic][i]
            # Aromatic oxygen (O1)
            mol.v[:WCLogP][i] = :O1
            continue
        end
        nbrs = collect(neighborkeys(mol, i))
        bonds = collect(neighboredgekeys(mol, i))
        if mol.v[:Degree][i] == 1
            # Hydroxyl (O2,5-12)
            nbr = nbrs[1]
            bond = bonds[1]
            if mol.v[:WCLogP][nbr] == :C5 && mol.v[:Charge][i] < 0
                # Acid (O12)
                mol.v[:WCLogP][i] = :O12
            elseif mol.v[:BondOrder][bond] == 2 || mol.v[:Charge][i] < 0
                # Carbonyl or oxide (O5-11)
                if mol.v[:WCLogP][nbr] == :C25
                    mol.v[:WCLogP][i] = :O8 # Aromatic carbonyl
                elseif mol.v[:WCLogP][nbr] == :C5
                    # Aliphatic carbonyl
                    cnbrs = neighborkeys(mol, nbr)
                    pop!(cnbrs, i)
                    cnbrs = collect(cnbrs)
                    if all(mol.v[:Symbol][cnbrs] .!= :C)
                        mol.v[:WCLogP][i] = :O11 # Carbonyl heteroatom
                    elseif any(mol.v[:Aromatic][cnbrs] .== true)
                        mol.v[:WCLogP][i] = :O10 # Carbonyl aromatic
                    else
                        mol.v[:WCLogP][i] = :O9 # Carbonyl aliphatic
                    end
                elseif mol.v[:Symbol][nbr] in (:O, :N)
                    mol.v[:WCLogP][i] = :O5 # O2? or N-oxide
                elseif mol.v[:Symbol][nbr] == :S
                    mol.v[:WCLogP][i] = :O6 # S-oxide
                else
                    mol.v[:WCLogP][i] = :O7 # Other oxide
                end
            else
                # Alcohol (O2)
                mol.v[:WCLogP][i] = :O2
            end
        elseif mol.v[:Degree][i] == 2
            # Ether (O3,4)
            if all(mol.v[:Aromatic][nbrs] .== false)
                mol.v[:WCLogP][i] = :O3 # Aliphatic ether
            else
                mol.v[:WCLogP][i] = :O4 # Aromatic ether
            end
        end
    end

    # Others
    for i in findall(mol.v[:WCLogP] .== :undef)
        if mol.v[:Symbol][i] in (:F, :Cl, :Br, :I)
            if mol.v[:Charge][i] == 0
                mol.v[:WCLogP][i] = mol.v[:Symbol][i]
            else
                mol.v[:WCLogP][i] = :Hal
            end
        elseif mol.v[:Symbol][i] == :P
            mol.v[:WCLogP][i] = :P
        elseif mol.v[:Symbol][i] == :S
            if mol.v[:Aromatic][i]
                mol.v[:WCLogP][i] = :S3
            elseif mol.v[:Charge][i] == 0
                mol.v[:WCLogP][i] = :S1
            else
                mol.v[:WCLogP][i] = :S2
            end
        elseif mol.v[:Symbol][i] == :H
            nbr = pop!(neighborkeys(mol, i))
            mol.v[:WCLogP][i] = wclogphydrogentype(mol, nbr)
        elseif mol.v[:Symbol][i] in P_BLOCK
            mol.v[:WCLogP][i] = :Me1
        elseif mol.v[:Symbol][i] in D_BLOCK
            mol.v[:WCLogP][i] = :Me2
        end
    end
end


function wclogphydrogentype(mol, i)
    if mol.v[:Symbol][i] == :C
        return :H1 # Hydrocarbon
    elseif mol.v[:Symbol][i] == :N
        return :H3 # Amine
    elseif mol.v[:Symbol][i] == :O
        nbr = pop!(neighborkeys(mol, i))
        if mol.v[:Symbol][nbr] == :N
            return :H3 # Hydroxyamine
        elseif mol.v[:Symbol][nbr] in (:O, :S)
            return :H4 # Peroxide, sulfoxide
        elseif mol.v[:Symbol][nbr] == :C
            if mol.v[:Pi][nbr] == 1 && !mol.v[:Aromatic][nbr]
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


function wclogpcontrib!(mol)
    mol.v[:WCLogPContrib] = zeros(atomcount(mol))
    for i in nodekeys(mol)
        hcnt = mol.v[:H_Count][i]
        contrib = get(WCLOGP_TABLE, string(mol.v[:WCLogP][i]), 0)
        if mol.v[:Symbol][i] == :H
            mol.v[:WCLogPContrib][i] = 0  # avoid double count
        elseif hcnt == 0
            mol.v[:WCLogPContrib][i] = contrib
        else
            htype = wclogphydrogentype(mol, i)
            hcontrib = WCLOGP_TABLE[string(htype)] * hcnt
            mol.v[:WCLogPContrib][i] = contrib + hcontrib
        end
    end
    return
end
