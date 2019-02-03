#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    trivialhydrogens,
    allhydrogens,
    largestcomponent,
    neutralize_acids!,
    neutralize_oniums!,
    depolarize!,
    triplebond_anion!,
    canonicalize!

# TODO: large conjugated system
# TODO: salts and waters should detected by functional group analysis
# TODO: Phosphate, diphosphate, sulfate, nitrate, acetate,
# maleate, fumarate, succinate, citrate, tartrate, oxalate,
# mesylate, tosylate, besylate,
# benzoate, gluconate

# TODO: hydrogensexplicit!(mol::MolGraph)


"""
    trivialhydrogens(mol::MolGraph) -> Set{Int}

Return a set of trivial hydrogen nodes (light hydrogen which is uncharged,
non-radical, non-stereospecific and attached to organic heavy atoms)
"""
function trivialhydrogens(mol::MolGraph)
    hs = Set{Int}()
    organicheavy = (
        :B, :C, :N, :O, :F, :Si, :P, :S, :Cl, :As, :Se, :Br, :I)
    for (i, a) in nodesiter(mol)
        if (a.symbol != :H || a.charge != 0 || a.multiplicity != 1
                || a.mass !== nothing)
            continue
        elseif a isa SmilesAtom && a.stereo !== nothing
            continue
        end
        nbrs = neighbors(mol, i)
        if length(nbrs) != 1
            continue
        end
        (nbr, bond) = pop!(nbrs)
        if bond isa SDFileBond && (bond.order != 1 || bond.notation != 0)
            continue
        elseif !in(getnode(mol, nbr).symbol, organicheavy)
            continue
        end
        push!(hs, i)
    end
    return hs
end


"""
    allhydrogens(mol::MolGraph) -> Set{Int}

Return a set of hydrogen nodes.
"""
function allhydrogens(mol::MolGraph)
    hs = Set{Int}()
    for (i, a) in nodesiter(mol)
        if a.symbol == :H
            push!(hs, i)
        end
    end
    return hs
end


"""
    largestcomponent(mol::MolGraph) -> Set{Int}

Return a set of nodes in the largest connected component.
"""
function largestcomponent(mol::MolGraph)
    # TODO: better way like python's max(iter, key=cmp)
    conn = connected_components(mol)
    sizemax = map(length, conn)
    largest = conn[argmax(sizemax)]
    return largest
end


"""
    neutralize_acids!(mol::VectorMol)

Neutralize oxo(thio) acids.

Note that this function edits `Atom` object fields directly. The molecular
property vector needs recalculation to apply the changes.
see [`canonicalize!`](@ref).
"""
function neutralize_acids!(mol::VectorMol)
    for o in findall((mol[:Symbol] .== :O)
            .* (mol[:Charge] .== -1) .* (mol[:Connectivity] .== 1))
        nbr = pop!(neighborkeys(mol, o))
        if mol[:Pi][nbr] == 1
            cnbrs = neighborkeys(mol, nbr)
            pop!(cnbrs, o)
            for cn in cnbrs
                if (mol[:Symbol][cn] in (:O, :S)
                        && mol[:Pi][cn] == 1 && mol[:Connectivity][cn] == 1)
                    oatom = getnode(mol, o)
                    oatom.charge = 0
                    break
                end
            end
        end
    end
end


"""
    neutralize_oniums!(mol::VectorMol)

Neutralize 1-3° oniums. Permanently charged quart-oniums are not neutralized.

Note that this function edits `Atom` object fields directly. The molecular
property vector needs recalculation to apply the changes.
see [`canonicalize!`](@ref).
"""
function neutralize_oniums!(mol::VectorMol)
    for o in findall((mol[:Charge] .== 1) .* (mol[:H_Count] .> 0))
        oatom = getnode(mol, o)
        oatom.charge = 0
    end
end


"""
    depolarize!(mol::VectorMol)

Depolarize oxo groups except in the case that polarization is required for
aromaticity.

Note that this function edits `Atom` object fields directly. The molecular
property vector needs recalculation to apply the changes.
see [`canonicalize!`](@ref).
"""
function depolarize!(mol::VectorMol)
    for o in findall((mol[:Symbol] .== :O) .* (mol[:Charge] .== -1))
        nbrs = neighborkeys(mol, o)
        @assert length(nbrs) == 1 "unexpected oxygen degree $(length(nbrs))"
        nbr = pop!(nbrs)
        if mol[:Charge][nbr] == 1 && !mol[:Aromatic][nbr]
            oatom = getnode(mol, o)
            oatom.charge = 0
            natom = getnode(mol, nbr)
            natom.charge = 0
        end
    end
end


"""
    triplebond_anion!(mol::VectorMol)

Canonicalize anions next to triple bonds (ex. [C-][N+]#N -> C=[N+]=[N-]).

Note that this function edits `Atom` object fields directly. The molecular
property vector needs recalculation to apply the changes.
see [`canonicalize!`](@ref).
"""
function triplebond_anion!(mol::VectorMol)
    # TODO: better function name
    for tb in findall(mol[:BondOrder] .== 3)
        tbond = getbond(mol, tb)
        for (f, s) in ((tbond.u, tbond.v), (tbond.v, tbond.u))
            nbrs = neighbors(mol, f)
            pop!(nbrs, s)
            if length(nbrs) != 1
                continue
            end
            (nbr, nbond) = pop!(nbrs)
            if mol[:Charge][nbr] == -1
                natom = getnode(mol, nbr)
                natom.charge = 0
                nbond.order = 2
                satom = getnode(mol, s)
                satom.charge = -1
                tbond.order = 2
            end
        end
    end
end


"""
    canonicalize!(mol::VectorMol)

Canonicalize molecule notation and apply the changes to the molecular property
vector.

- Neutralize oxo acid, 1-3° ammonium and polarized carbonyls except in the
  case that polarization is required for aromaticity.
- Canonicalize anions next to triple bonds (ex. [C-][N+]#N -> C=[N+]=[N-])
"""
function canonicalize!(mol::VectorMol)
    aromatic!(mol)
    neutralizeacids!(mol)
    neutralizeoniums!(mol)
    depolarize!(mol)
    triplebondanion!(mol)
    aromatic!(mol, recalculate=true)
end
