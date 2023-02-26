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
    wclogptype(mol::MolGraph)

Return Wildman-Crippen LogP atom types.
"""
function wclogptype(mol::MolGraph)
    atomtypes = init_node_descriptor(Symbol, mol)
    fill!(atomtypes, :undef)
    atomsymbol_ = atom_symbol(mol)
    charge_ = charge(mol)
    bondorder_ = bond_order(mol)
    heavyatoms_ = heavy_atoms(mol)
    hydrogens_ = total_hydrogens(mol)
    pielectron_ = pi_electron(mol)
    hybridization_ = hybridization(mol)
    isaromatic_ = is_aromatic(mol)
    isaromaticbond_ = is_edge_aromatic(mol)

    for i in vertices(mol)
        # Carbons
        if atomsymbol_[i] === :C && hybridization_[i] === :sp3
            # Aliphatic sp3 hybrid (C1-4,8-12,27)
            nbrs = neighbors(mol, i)
            if !isempty(setdiff(atomsymbol_[nbrs], ALIPH_HETERO))
                atomtypes[i] = :C27  # Adjacent to inorganic atoms
            elseif all(isaromatic_[nbrs] .=== false)
                if all(atomsymbol_[nbrs] .=== :C)
                    # Aliphatic carbon (C1,2)
                    atomtypes[i] = heavyatoms_[i] < 3 ? :C1 : :C2
                else
                    # Adjacent to heteroatoms (C3,4)
                    atomtypes[i] = heavyatoms_[i] < 3 ? :C3 : :C4
                end
            else
                # Adjacent to aromatic atoms (C8-12)
                if heavyatoms_[i] == 1
                    atomtypes[i] = atomsymbol_[nbrs[1]] === :C ? :C8 : :C9
                elseif heavyatoms_[i] == 2
                    atomtypes[i] = :C10
                elseif heavyatoms_[i] == 3
                    atomtypes[i] = :C11
                elseif heavyatoms_[i] == 4
                    atomtypes[i] = :C12
                end
            end
        elseif atomsymbol_[i] === :C && isaromatic_[i]
            # Aromatic (C13-25)
            if hydrogens_[i] == 1
                atomtypes[i] = :C18  # Aromatic non-substituted
                continue
            end
            subst = -1
            substbond = -1
            aromcnt = 0
            for nbr in neighbors(mol, i)
                if isaromaticbond_[edge_rank(mol, i, nbr)]
                    aromcnt += 1
                else
                    subst = nbr
                    substbond = edge_rank(mol, i,nbr)
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
        elseif atomsymbol_[i] === :C && hybridization_[i] === :sp2
            # Aliphatic sp2 hybrid (C5-6,26)
            arom = 0
            het = 0
            for nbr in neighbors(mol, i)
                if bondorder_[edge_rank(mol, i, nbr)] == 2 && atomsymbol_[nbr] !== :C
                    het += 1
                elseif isaromatic_[nbr]
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
        elseif atomsymbol_[i] === :C && hybridization_[i] === :sp
            # Aliphatic sp hybrid (C6-7)
            bonds = [edge_rank(mol, i, nbr) for nbr in neighbors(mol, i)]
            if any(bondorder_[bonds] .== 3)
                atomtypes[i] = :C7 # Alkyne, Nitrile
            else
                atomtypes[i] = :C6 # Allene
            end

        # Nitrogens
        elseif atomsymbol_[i] === :N && hybridization_[i] === :sp3
            if charge_[i] > 0
                # Ammonium (N10, N13)
                if heavyatoms_[i] == 4
                    atomtypes[i] = :N13  # Quart-ammonium
                elseif heavyatoms_[i] < 4
                    atomtypes[i] = :N10  # Protonated amine
                end
            elseif charge_[i] == 0
                # Aliphatic amine (N1,2,7)
                if heavyatoms_[i] == 1
                    atomtypes[i] = :N1
                elseif heavyatoms_[i] == 2
                    atomtypes[i] = :N2
                elseif heavyatoms_[i] == 3
                    atomtypes[i] = :N7
                end
            end
        elseif atomsymbol_[i] === :N && isaromatic_[i]
            # Aromatic amine (N11,12)
            if charge_[i] > 0
                atomtypes[i] = :N12 # Protonated aromatic
            elseif charge_[i] == 0
                atomtypes[i] = :N11 # Unprotonated aromatic
            end
        elseif atomsymbol_[i] === :N && hybridization_[i] === :sp2
            if pielectron_[i] == 2
                nbrs = neighbors(mol, i)
                if any(isaromatic_[nbrs])
                    # Amine adjacent to aromatic (N3,4,8)
                    if heavyatoms_[i] == 1
                        atomtypes[i] = :N3
                    elseif heavyatoms_[i] == 2
                        atomtypes[i] = :N4
                    elseif heavyatoms_[i] == 3
                        atomtypes[i] = :N8
                    end
                else
                    # Aliphatic amine (N1,2,7)
                    if heavyatoms_[i] == 1
                        atomtypes[i] = :N1
                    elseif heavyatoms_[i] == 2
                        atomtypes[i] = :N2
                    elseif heavyatoms_[i] == 3
                        atomtypes[i] = :N7
                    end
                end
            elseif pielectron_[i] == 1
                # Imine (N5,6)
                if charge_[i] > 0
                    atomtypes[i] = :N13  # Quart-ammonium
                elseif charge_[i] == 0
                    if heavyatoms_[i] == 1
                        atomtypes[i] = :N5
                    elseif heavyatoms_[i] == 2
                        atomtypes[i] = :N6
                    end
                end
            end
        elseif atomsymbol_[i] === :N && hybridization_[i] === :sp
            # sp nitrogen (N9,13,14)
            adjatoms = [atomsymbol_[nbr] for nbr in neighbors(mol, i)]
            if charge_[i] == 1 && issetequal(adjatoms, [:N, :N])
                atomtypes[i] = :N14  # Other ionized nitrogen (Azide)
            elseif charge_[i] == 1 && issetequal(adjatoms, [:C, :N])
                atomtypes[i] = :N13  # Quart-ammonium (Diazo)
            elseif charge_[i] == 0 && :C in adjatoms
                atomtypes[i] = :N9 # Nitrile
            end

        # Oxygens
        elseif atomsymbol_[i] === :O && hybridization_[i] === :sp3
            if hydrogens_[i] > 0
                atomtypes[i] = :O2  # Alcohol (O2)
            elseif charge_[i] == 0 && heavyatoms_[i] == 2
                atomtypes[i] = :O3  # Aliphatic ether
            elseif charge_[i] < 0 && heavyatoms_[i] == 1
                # Oxide (O5,6,7,12)
                nbr = neighbors(mol, i)[1]
                if atomsymbol_[nbr] in (:O, :N)
                    atomtypes[i] = :O5 # O2 or N-oxide
                elseif atomsymbol_[nbr] == :S
                    atomtypes[i] = :O6 # S-oxide
                elseif atomsymbol_[nbr] == :C
                    cnbrsyms = [atomsymbol_[n] for n in neighbors(mol, nbr)]
                    if sum(cnbrsyms .=== :O) == 2
                        atomtypes[i] = :O12 # Acid (O12)
                    end
                end
                if atomtypes[i] === :undef
                    atomtypes[i] = :O7 # Other oxide
                end
            end
        elseif atomsymbol_[i] === :O && isaromatic_[i]
            atomtypes[i] = :O1  # Aromatic oxygen (O1)
        elseif atomsymbol_[i] === :O && hybridization_[i] === :sp2 && charge_[i] == 0
            if heavyatoms_[i] == 2
                if any(isaromatic_[neighbors(mol, i)])
                    atomtypes[i] = :O4  # Aromatic ether
                else
                    atomtypes[i] = :O3  # Aliphatic ether
                end
            elseif pielectron_[i] == 2
                atomtypes[i] = :O2  # Alcohol
            else
                # Carbonyl
                c = neighbors(mol, i)[1]
                _, cnbrs = edge_neighbors(mol, i, c)
                if isaromatic_[c]
                    atomtypes[i] = :O8  # Aromatic carbonyl
                elseif all(atomsymbol_[cnbrs] .!== :C)
                    atomtypes[i] = :O11  # Carbonyl heteroatom
                elseif any(isaromatic_[cnbrs])
                    atomtypes[i] = :O10  # Carbonyl aromatic
                else
                    atomtypes[i] = :O9  # Carbonyl aliphatic
                end
            end

        # Others
        elseif atomsymbol_[i] in (:F, :Cl, :Br, :I)  # ?  && atomtypes[i] === :undef
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
            # Keep undef (defined by adjacent wclogphydrogentype)
            continue
        elseif atomsymbol_[i] in P_BLOCK
            atomtypes[i] = :Me1
        elseif atomsymbol_[i] in D_BLOCK
            atomtypes[i] = :Me2
        end

        if atomtypes[i] === :undef
            if atomsymbol_[i] === :C
                atomtypes[i] = :CS  # Other carbons
            elseif atomsymbol_[i] === :O
                atomtypes[i] = :OS  # Other oxygens
            elseif atomsymbol_[i] === :N
                if charge_[i] < 0
                    atomtypes[i] = :N14  # Other ionized nitrogen (Nitride)
                else
                    atomtypes[i] = :NS  # Other nitrogens
                end
            end
        end
    end
    return atomtypes
end


function wclogphydrogentype(mol::MolGraph)
    atomtypes = init_node_descriptor(Symbol, mol)
    fill!(atomtypes, :undef)
    atomsymbol_ = atom_symbol(mol)
    degree_ = degree(mol)
    heavyatoms_ = heavy_atoms(mol)
    hydrogens_ = total_hydrogens(mol)
    hybridization_ = hybridization(mol)
    isaromatic_ = is_aromatic(mol)

    for i in vertices(mol)
        hydrogens_[i] == 0 && continue
        if atomsymbol_[i] === :H && degree_[i] == 1
            atomtypes[i] = :HS  # Proton or hydride (But not meaningful)
        elseif atomsymbol_[i] in (:H, :C)
            atomtypes[i] = :H1  # Hydrocarbon or molecular hydrogen
        elseif atomsymbol_[i] === :N
            atomtypes[i] = :H3  # Amine
        elseif atomsymbol_[i] === :O
            if heavyatoms_[i] == 0
                atomtypes[i] = :H2  # Alcohol (H2O)
                continue
            end
            heavy = filter(x -> atomsymbol_[x] !== :H, neighbors(mol, i))[1]
            if atomsymbol_[heavy] === :N
                atomtypes[i] = :H3  # Hydroxyamine
            elseif atomsymbol_[heavy] in (:O, :S)
                atomtypes[i] = :H4  # Peroxide, sulfoxide
            elseif atomsymbol_[heavy] === :C
                if hybridization_[heavy] === :sp2 && !isaromatic_[heavy]
                    atomtypes[i] = :H4  # Acid
                else
                    atomtypes[i] = :H2  # Alcohol
                end
            else
                atomtypes[i] = :H2  # Alcohol (Other hydroxyl)
            end
        else
            atomtypes[i] = :H2  # Other hydrate
        end
    end
    return atomtypes
end


function wclogpcontrib(mol::MolGraph)
    contrib = init_node_descriptor(Float64, mol)
    atomsymbol_ = atom_symbol(mol)
    heavyatoms_ = heavy_atoms(mol)
    hydrogens_ = total_hydrogens(mol)
    logptypes_ = wclogptype(mol)
    htypes_ = wclogphydrogentype(mol)

    for i in vertices(mol)
        cont = get(WCLOGP_TABLE, string(logptypes_[i]), 0.0)
        hcont = get(WCLOGP_TABLE, string(htypes_[i]), 0.0)
        if atomsymbol_[i] === :H
            if heavyatoms_[i] > 0
                contrib[i] = 0.0  # avoid double count
            else
                contrib[i] = hcont
            end
        else
            contrib[i] = cont + hcont * hydrogens_[i]
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
wclogp(mol::MolGraph) = reduce(+, wclogpcontrib(mol); init=0)
wclogp(mol::MolGraph, digits) = round(wclogp(mol), digits=digits)