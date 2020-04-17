#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    wclogptype, wclogphydrogentype,
    wclogpcontrib, wclogp


const WCLOGP_TABLE = YAML.load(open(
    joinpath(dirname(@__FILE__), "..", "assets", "const", "wclogp.yaml")
))["logP"]

const ALIPH_HETERO = Set([:H, :C, :N, :O, :P, :S, :F, :Cl, :Br, :I])
const AROM_HETERO = Dict(
    :H  => :C18, :C => :C21, :N => :C22, :O => :C23, :S => :C24,
    :F => :C14, :Cl => :C15, :Br => :C16, :I => :C17
)
# TODO:
const P_BLOCK = (:Al, :B, :Si, :Ga, :Ge, :As, :Se, :Sn, :Te, :Pb,
                 :Ne, :Ar, :Kr, :Xe, :Rn)
const D_BLOCK = (:Fe, :Co, :Cu, :Zn, :Tc, :Cd, :Pt, :Au, :Hg, :Gd)


"""
    wclogptype(mol::GraphMol)

Return Wildman-Crippen LogP atom types.
"""
@cache function wclogptype(mol::GraphMol)
    atomtypes = fill(:undef, nodecount(mol))
    atomsymbol_ = atomsymbol(mol)
    charge_ = charge(mol)
    bondorder_ = bondorder(mol)
    heavyatomcount_ = heavyatomcount(mol)
    hcount_ = hcount(mol)
    pielectron_ = pielectron(mol)
    hybridization_ = hybridization(mol)
    isaromatic_ = isaromatic(mol)
    isaromaticbond_ = isaromaticbond(mol)
    # Carbons
    for i in findall(atomsymbol_ .=== :C)
        if hybridization_[i] === :sp3
            # Aliphatic sp3 hybrid (C1-4,8-12,27)
            nbrs = collect(adjacencies(mol, i))
            if !isempty(setdiff(atomsymbol_[nbrs], ALIPH_HETERO))
                atomtypes[i] = :C27  # Adjacent to inorganic atoms
            elseif all(isaromatic_[nbrs] .=== false)
                if all(atomsymbol_[nbrs] .=== :C)
                    # Aliphatic carbon (C1,2)
                    atomtypes[i] = heavyatomcount_[i] < 3 ? :C1 : :C2
                else
                    # Adjacent to heteroatoms (C3,4)
                    atomtypes[i] = heavyatomcount_[i] < 3 ? :C3 : :C4
                end
            else
                # Adjacent to aromatic atoms (C8-12)
                if heavyatomcount_[i] == 1
                    atomtypes[i] = atomsymbol_[nbrs[1]] === :C ? :C8 : :C9
                elseif heavyatomcount_[i] == 2
                    atomtypes[i] = :C10
                elseif heavyatomcount_[i] == 3
                    atomtypes[i] = :C11
                elseif heavyatomcount_[i] == 4
                    atomtypes[i] = :C12
                end
            end
        elseif isaromatic_[i]
            # Aromatic (C13-25)
            if hcount_[i] == 1
                atomtypes[i] = :C18  # Aromatic non-substituted
                continue
            end
            subst = -1
            substbond = -1
            aromcnt = 0
            for (inc, adj) in neighbors(mol, i)
                if isaromaticbond_[inc]
                    aromcnt += 1
                else
                    subst = adj
                    substbond = inc
                end
            end
            if aromcnt == 3
                atomtypes[i] = :C19 # Bridgehead
            elseif !haskey(AROM_HETERO, atomsymbol_[subst])
                atomtypes[i] = :C13 # Inorganic substituent
            elseif isaromatic_[subst]
                atomtypes[i] = :C20 # Aromatic substituent
            elseif bondorder_[substbond] == 2
                atomtypes[i] = :C25 # Double bond substituent
            else
                # Typical substituent (C14-17,21-24)
                atomtypes[i] = AROM_HETERO[atomsymbol_[subst]]
            end
        elseif hybridization_[i] === :sp2
            # Aliphatic sp2 hybrid (C5-6,26)
            arom = 0
            het = 0
            for (inc, adj) in neighbors(mol, i)
                if bondorder_[inc] == 2 && atomsymbol_[adj] !== :C
                    het += 1
                elseif isaromatic_[adj]
                    arom += 1
                end
            end
            if het > 0
                atomtypes[i] = :C5 # Double bond to non-C
            elseif arom > 0
                atomtypes[i] = :C26 # Adjacent to aromatic
            else
                atomtypes[i] = :C6 # Double bond to C
            end
        elseif hybridization_[i] === :sp
            # Aliphatic sp hybrid (C6-7)
            bonds = collect(incidences(mol, i))
            if any(bondorder_[bonds] .== 3)
                atomtypes[i] = :C7 # Alkyne, Nitrile
            else
                atomtypes[i] = :C6 # Allene
            end
        end
        if atomtypes[i] === :undef
            atomtypes[i] = :CS # Not found
        end
    end

    # Nitrogens
    for i in findall(atomsymbol_ .=== :N)
        if hybridization_[i] === :sp3
            if charge_[i] > 0
                # Ammonium (N10, N13)
                if heavyatomcount_[i] == 4
                    atomtypes[i] = :N13  # Quart-ammonium
                elseif heavyatomcount_[i] < 4
                    atomtypes[i] = :N10  # Protonated amine
                end
            elseif charge_[i] == 0
                # Aliphatic amine (N1,2,7)
                if heavyatomcount_[i] == 1
                    atomtypes[i] = :N1
                elseif heavyatomcount_[i] == 2
                    atomtypes[i] = :N2
                elseif heavyatomcount_[i] == 3
                    atomtypes[i] = :N7
                end
            end
        elseif isaromatic_[i]
            # Aromatic amine (N11,12)
            if charge_[i] > 0
                atomtypes[i] = :N12 # Protonated aromatic
            elseif charge_[i] == 0
                atomtypes[i] = :N11 # Unprotonated aromatic
            end
        elseif hybridization_[i] === :sp2
            if pielectron_[i] == 2
                adjs = collect(adjacencies(mol, i))
                if any(isaromatic_[adjs])
                    # Amine adjacent to aromatic (N3,4,8)
                    if heavyatomcount_[i] == 1
                        atomtypes[i] = :N3
                    elseif heavyatomcount_[i] == 2
                        atomtypes[i] = :N4
                    elseif heavyatomcount_[i] == 3
                        atomtypes[i] = :N8
                    end
                else
                    # Aliphatic amine (N1,2,7)
                    if heavyatomcount_[i] == 1
                        atomtypes[i] = :N1
                    elseif heavyatomcount_[i] == 2
                        atomtypes[i] = :N2
                    elseif heavyatomcount_[i] == 3
                        atomtypes[i] = :N7
                    end
                end
            elseif pielectron_[i] == 1
                # Imine (N5,6)
                if charge_[i] > 0
                    atomtypes[i] = :N13  # Quart-ammonium
                elseif charge_[i] == 0
                    if heavyatomcount_[i] == 1
                        atomtypes[i] = :N5
                    elseif heavyatomcount_[i] == 2
                        atomtypes[i] = :N6
                    end
                end
            end
        elseif hybridization_[i] === :sp
            # sp nitrogen (N9,13,14)
            adjs = adjacencies(mol, i)
            adjatoms = [atomsymbol_[adj] for adj in adjs]
            if charge_[i] == 1 && issetequal(adjatoms, [:N, :N])
                atomtypes[i] = :N14  # Other ionized nitrogen (Azide)
            elseif charge_[i] == 1 && issetequal(adjatoms, [:C, :N])
                atomtypes[i] = :N13  # Quart-ammonium (Diazo)
            elseif charge_[i] == 0 && :C in adjatoms
                atomtypes[i] = :N9 # Nitrile
            end
        end
        if atomtypes[i] === :undef
            if charge_[i] < 0
                atomtypes[i] = :N14  # Other ionized nitrogen (Nitride)
            else
                atomtypes[i] = :NS  # Others
            end
        end
    end

    # Oxygens
    for i in findall(atomsymbol_ .=== :O)
        if hybridization_[i] === :sp3
            if hcount_[i] > 0
                atomtypes[i] = :O2  # Alcohol (O2)
            elseif charge_[i] == 0 && heavyatomcount_[i] == 2
                atomtypes[i] = :O3  # Aliphatic ether
            elseif charge_[i] < 0 && heavyatomcount_[i] == 1
                # Oxide (O5,6,7,12)
                adj = pop!(adjacencies(mol, i))
                if atomsymbol_[adj] in (:O, :N)
                    atomtypes[i] = :O5 # O2 or N-oxide
                elseif atomsymbol_[adj] == :S
                    atomtypes[i] = :O6 # S-oxide
                elseif atomsymbol_[adj] == :C
                    cadjs = adjacencies(mol, adj)
                    cadjatoms = [atomsymbol_[cadj] for cadj in cadjs]
                    if sum(cadjatoms .=== :O) == 2
                        atomtypes[i] = :O12 # Acid (O12)
                    end
                end
                if atomtypes[i] === :undef
                    atomtypes[i] = :O7 # Other oxide
                end
            end
        elseif isaromatic_[i]
            atomtypes[i] = :O1  # Aromatic oxygen (O1)
        elseif hybridization_[i] === :sp2 && charge_[i] == 0
            if heavyatomcount_[i] == 2
                adjs = collect(adjacencies(mol, i))
                if any(isaromatic_[adjs])
                    atomtypes[i] = :O4  # Aromatic ether
                else
                    atomtypes[i] = :O3  # Aliphatic ether
                end
            elseif pielectron_[i] == 2
                atomtypes[i] = :O2  # Alcohol
            else
                # Carbonyl
                adj = pop!(adjacencies(mol, i))
                cadjs = adjacencies(mol, adj)
                pop!(cadjs, i)
                cadjs = collect(cadjs)
                if isaromatic_[adj]
                    atomtypes[i] = :O8  # Aromatic carbonyl
                elseif all(atomsymbol_[cadjs] .!== :C)
                    atomtypes[i] = :O11  # Carbonyl heteroatom
                elseif any(isaromatic_[cadjs])
                    atomtypes[i] = :O10  # Carbonyl aromatic
                else
                    atomtypes[i] = :O9  # Carbonyl aliphatic
                end
            end
        end
        if atomtypes[i] === :undef
            atomtypes[i] = :OS  # Others
        end
    end

    # Others
    for i in findall(atomtypes .=== :undef)
        if atomsymbol_[i] in (:F, :Cl, :Br, :I)
            if charge_[i] == 0
                atomtypes[i] = atomsymbol_[i]
            else
                atomtypes[i] = :Hal  # Ionic halogens
            end
        elseif atomsymbol_[i] === :P
            atomtypes[i] = :P
        elseif atomsymbol_[i] === :S
            if isaromatic_[i]
                atomtypes[i] = :S3  # Aromatic sulfur
            elseif charge_[i] === 0
                atomtypes[i] = :S1  # Aliphatic sulfur
            else
                atomtypes[i] = :S2  # Ionic sulfur
            end
        elseif atomsymbol_[i] === :H
            adjs = adjacencies(mol, i)
            if isempty(adjs)
                atomtypes[i] = :HS  # Proton or hydride (But not meaningful)
            else
                adj = pop!(adjs)
                atomtypes[i] = wclogphydrogentype(mol, adj)
            end
        elseif atomsymbol_[i] in P_BLOCK
            atomtypes[i] = :Me1
        elseif atomsymbol_[i] in D_BLOCK
            atomtypes[i] = :Me2
        end
    end
    return atomtypes
end


function wclogphydrogentype(mol::GraphMol, i)
    atomsymbol_ = atomsymbol(mol)
    hybridization_ = hybridization(mol)
    isaromatic_ = isaromatic(mol)
    if atomsymbol_[i] in (:H, :C)
        return :H1  # Hydrocarbon
    elseif atomsymbol_[i] === :N
        return :H3  # Amine
    elseif atomsymbol_[i] === :O
        adjs = adjacencies(mol, i)
        heavy = [adj for adj in adjs if atomsymbol_[adj] !== :H]
        isempty(heavy) && return :H2  # Alcohol (H2O)
        a = pop!(heavy)
        if atomsymbol_[a] === :N
            return :H3  # Hydroxyamine
        elseif atomsymbol_[a] in (:O, :S)
            return :H4  # Peroxide, sulfoxide
        elseif atomsymbol_[a] === :C
            if hybridization_[a] === :sp2 && !isaromatic_[a]
                return :H4  # Acid
            else
                return :H2  # Alcohol
            end
        else
            return :H2  # Alcohol (Other hydroxyl)
        end
    else
        return :H2  # Other hydrate
    end
end


function wclogpcontrib(mol::GraphMol)
    contrib = zeros(nodecount(mol))
    atomsymbol_ = atomsymbol(mol)
    hcount_ = hcount(mol)
    atomtypes = wclogptype(mol)
    for i in 1:nodecount(mol)
        cont = get(WCLOGP_TABLE, string(atomtypes[i]), 0)
        if atomsymbol_[i] == :H
            contrib[i] = 0  # avoid double count
        elseif hcount_[i] == 0
            contrib[i] = cont
        else
            htype = wclogphydrogentype(mol, i)
            hcont = WCLOGP_TABLE[string(htype)] * hcount_[i]
            contrib[i] = cont + hcont
        end
    end
    return contrib
end


"""
    wclogp(mol::GraphMol) -> Float64

Return predicted logP value calculated by using Wildman and Crippen method.

# Reference

1. Wildman, S. A. and Crippen, G. M. (1999). Prediction of Physicochemical
   Parameters by Atomic Contributions. Journal of Chemical Information and
   Modeling, 39(5), 868–873. https://doi.org/10.1021/ci990307l
"""
wclogp(mol::GraphMol; digits=2
    ) = round(reduce(+, wclogpcontrib(mol); init=0), digits=digits)
