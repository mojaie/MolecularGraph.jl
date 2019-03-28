#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    trivialhydrogens,
    all_hydrogens,
    make_hydrogens_implicit,
    make_hydrogens_explicit,
    largest_component_nodes,
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
    trivialhydrogens(mol::VectorMol) -> Set{Int}

Return a set of trivial hydrogen nodes (light hydrogens which are uncharged,
non-radical, non-stereospecific and attached to organic heavy atoms)
"""
function trivialhydrogens(mol::VectorMol)
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
    all_hydrogens(mol::VectorMol) -> Set{Int}

Return a set of hydrogen nodes.
"""
function all_hydrogens(mol::VectorMol)
    hs = Set{Int}()
    for (i, a) in nodesiter(mol)
        if a.symbol == :H
            push!(hs, i)
        end
    end
    return hs
end


"""
    make_hydrogens_implicit(mol::VectorMol) -> VectorMol

Return molecule whose hydrogen nodes are removed. If option `all` is set to
false, only trivial hydrogens are removed (see [`trivialhydrogens`](@ref)).
"""
function make_hydrogens_implicit(mol::VectorMol; all=true)
    hydrogens = all ? all_hydrogens : trivialhydrogens
    ns = setdiff(nodeset(mol), hydrogens(mol))
    return nodesubgraph(mol, ns)
end


"""
    make_hydrogens_explicit(mol::VectorMol) -> VectorMol

Return molecule whose hydrogens are fully attached. If option `all` is set to
false, only trivial hydrogens are removed (see [`trivialhydrogens`](@ref)).
"""
function make_hydrogens_explicit(mol::VectorMol)
    newmol = vectormol(mol)
    ncnt = nodecount(mol)
    impl = implicithcount(mol)
    for (n, node) in nodesiter(mol)
        for i in 1:impl[n]
            ncnt += 1
            updatenode!(newmol, nodetype(mol)(:H), ncnt)
            updateedge!(newmol, edgetype(mol)(n, ncnt))
        end
    end
    return newmol
end


"""
    largest_component_nodes(mol::VectorMol) -> Set{Int}

Return a set of nodes in the largest connected component.
"""
function largest_component_nodes(mol::VectorMol)
    # TODO: better way like python's max(iter, key=cmp)
    conn = connected_components(mol)
    sizemax = map(length, conn)
    largest = conn[argmax(sizemax)]
    return largest
end


"""
    largestcomponent(mol::VectorMol) -> VectorMol

Return largest connected component of the molecular graph.
"""
function largestcomponent(mol::VectorMol)
    ns = largest_component_nodes(mol)
    return nodesubgraph(mol, ns)
end


"""
    neutralize_acids!(mol::VectorMol)

Neutralize oxo(thio) acids.

Note that this function edits `Atom` object fields directly. The molecular
property vector needs recalculation to apply the changes.
see [`canonicalize!`](@ref).
"""
function neutralize_acids!(mol::VectorMol)
    symbols = atomsymbol(mol)
    charges = charge(mol)
    conn = connectivity(mol)
    pies = pielectron(mol)
    for o in findall((symbols .== :O) .* (charges .== -1) .* (conn .== 1))
        nbr = pop!(adjacencies(mol, o))
        if pies[nbr] == 1
            cnbrs = adjacencies(mol, nbr)
            pop!(cnbrs, o)
            for cn in cnbrs
                if symbols[cn] in (:O, :S) && pies[cn] == 1 && conn[cn] == 1
                    oatom = setcharge(getnode(mol, o), 0)
                    updatenode!(mol, oatom, o)
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
    for o in findall((mol[:charge] .== 1) .* (mol[:hcount] .> 0))
        oatom = setcharge(getnode(mol, o), 0)
        updatenode!(mol, oatom, o)
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
    symbols = atomsymbol(mol)
    charges = charge(mol)
    isarom = isaromatic(mol)
    for o in findall((symbols .== :O) .* (charges .== -1))
        nbrs = collect(neighbors(mol, o))
        @assert length(nbrs) == 1 "unexpected oxygen degree $(length(nbrs))"
        (nbr, b) = pop!(nbrs)
        if charges[nbr] == 1 && !isarom[nbr]
            oatom = setcharge(getnode(mol, o), 0)
            updatenode!(mol, oatom, o)
            natom = setcharge(getnode(mol, nbr), 0)
            updatenode!(mol, natom, nbr)
            bond = setorder(getedge(mol, b), 2)
            updateedge!(mol, bond, b)
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
    for tb in findall(mol[:bondorder] .== 3)
        tbond = getbond(mol, tb)
        for (f, s) in ((tbond.u, tbond.v), (tbond.v, tbond.u))
            nbrs = copy(neighbors(mol, f))
            pop!(nbrs, s)
            if length(nbrs) != 1
                continue
            end
            (nbr, nb) = pop!(nbrs)
            if mol[:charge][nbr] == -1
                natom = setcharge(getnode(mol, nbr), 0)
                updatenode!(mol, natom, nbr)
                satom = setcharge(getnode(mol, s), -1)
                updatenode!(mol, satom, s)
                tbond = setorder(getedge(mol, tb), 2)
                updateedge!(mol, tbond, tb)
                nbond = setorder(getedge(mol, nb), 2)
                updateedge!(mol, nbond, nb)
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
    neutralizeacids!(mol)
    neutralizeoniums!(mol)
    depolarize!(mol)
    triplebondanion!(mol)
    return
end

function canonicalize(mol::VectorMol)
    newmol = vectormol(mol)
    canonicalize!(newmol)
    return newmol
end
