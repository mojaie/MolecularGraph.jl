#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    trivialhydrogens,
    all_hydrogens,
    make_hydrogens_implicit!,
    make_hydrogens_implicit,
    make_hydrogens_explicit!,
    make_hydrogens_explicit,
    largest_component_nodes,
    largestcomponent!,
    largestcomponent,
    neutralize_acids!,
    neutralize_oniums!,
    depolarize!,
    triplebond_anion!,
    canonicalize!,
    canonicalize

# TODO: large conjugated system
# TODO: salts and waters should detected by functional group analysis
# TODO: Phosphate, diphosphate, sulfate, nitrate, acetate,
# maleate, fumarate, succinate, citrate, tartrate, oxalate,
# mesylate, tosylate, besylate,
# benzoate, gluconate


"""
    trivialhydrogens(mol::MolGraph) -> Set{Int}

Return a set of trivial hydrogen nodes (light hydrogens which are uncharged,
non-radical, non-stereospecific and attached to organic heavy atoms)
"""
function trivialhydrogens(mol::MolGraph)
    hs = Set{Int}()
    organic_heavy = (
        :B, :C, :N, :O, :F, :Si, :P, :S, :Cl, :As, :Se, :Br, :I)
    for (i, a) in nodesiter(mol)
        if (a.symbol != :H || a.charge != 0 || a.multiplicity != 1
                || a.mass !== nothing)
            continue
        elseif a isa SmilesAtom && a.stereo !== nothing
            continue
        end
        nbrs = collect(neighbors(mol, i))
        if length(nbrs) != 1
            continue
        end
        (nbr, bond) = pop!(nbrs)
        if bond isa SDFileBond && (bond.order != 1 || bond.notation != 0)
            continue
        elseif !in(getnode(mol, nbr).symbol, organic_heavy)
            continue
        end
        push!(hs, i)
    end
    return hs
end


"""
    all_hydrogens(mol::MolGraph) -> Set{Int}

Return a set of hydrogen nodes.
"""
function all_hydrogens(mol::MolGraph)
    hs = Set{Int}()
    for (i, a) in nodesiter(mol)
        if a.symbol == :H
            push!(hs, i)
        end
    end
    return hs
end


"""
    make_hydrogens_implicit!(mol::MolGraph) -> GeneralMapMol

Return molecule whose hydrogen nodes are removed. If option `all` is set to
false, only trivial hydrogens are removed (see [`trivialhydrogens`](@ref)).
"""
function make_hydrogens_implicit!(mol::MolGraph; all=true)
    hydrogens = all ? all_hydrogens : trivialhydrogens
    ns = setdiff(nodeset(mol), hydrogens(mol))
    return newgraph(nodesubgraph(mol, ns))
end

make_hydrogens_implicit(m::MolGraph; kwargs...) = make_hydrogens_implicit!(
    deepcopy(m); kwargs...)


"""
    make_hydrogens_explicit(mol::VectorMol) -> GeneralMapMol

Return molecule whose hydrogens are fully attached. If option `all` is set to
false, only trivial hydrogens are removed (see [`trivialhydrogens`](@ref)).
"""
function make_hydrogens_explicit!(mol::VectorMol)
    newmol = newgraph(mol)
    ncnt = nodecount(mol)
    ecnt = edgecount(mol)
    for (n, node) in nodesiter(mol)
        hcnt = mol[:Connectivity][n] - mol[:Degree][n]
        for i in 1:hcnt
            ncnt += 1
            ecnt += 1
            updatenode!(newmol, nodetype(mol)(:H), ncnt)
            updateedge!(newmol, edgetype(mol)(n, ncnt), ecnt)
        end
    end
    return newmol
end

make_hydrogens_explicit(m::VectorMol) = make_hydrogens_explicit!(deepcopy(m))


"""
    largest_component_nodes(mol::MolGraph) -> Set{Int}

Return a set of nodes in the largest connected component.
"""
function largest_component_nodes(mol::MolGraph)
    # TODO: better way like python's max(iter, key=cmp)
    conn = connected_components(mol)
    sizemax = map(length, conn)
    largest = conn[argmax(sizemax)]
    return largest
end


"""
    largestcomponent(mol::MolGraph) -> GeneralMapMol

Return largest connected component of the molecular graph.
"""
function largestcomponent!(mol::MolGraph)
    ns = largest_component_nodes(mol)
    return newgraph(nodesubgraph(mol, ns))
end

largestcomponent(mol::MolGraph) = largestcomponent!(deepcopy(mol))


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
        nbr = pop!(neighborset(mol, o))
        if mol[:Pi][nbr] == 1
            cnbrs = neighborset(mol, nbr)
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
        nbrs = collect(neighbors(mol, o))
        @assert length(nbrs) == 1 "unexpected oxygen degree $(length(nbrs))"
        (nbr, b) = pop!(nbrs)
        if mol[:Charge][nbr] == 1 && !mol[:Aromatic][nbr]
            oatom = getnode(mol, o)
            oatom.charge = 0
            natom = getnode(mol, nbr)
            natom.charge = 0
            bond = getedge(mol, b)
            bond.order = 2
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
            nbrs = copy(neighbors(mol, f))
            pop!(nbrs, s)
            if length(nbrs) != 1
                continue
            end
            (nbr, nb) = pop!(nbrs)
            if mol[:Charge][nbr] == -1
                natom = getnode(mol, nbr)
                natom.charge = 0
                satom = getnode(mol, s)
                satom.charge = -1
                nbond = getedge(mol, nb)
                nbond.order = 2
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
    return newgraph(mol)
end

canonicalize(mol::VectorMol) = canonicalize!(deepcopy(mol))
