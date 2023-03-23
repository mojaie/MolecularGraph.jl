#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    kekulize, kekulize!,
    removable_hydrogens, all_hydrogens,
    remove_hydrogens!, add_hydrogens!,
    largest_component_nodes, extract_largest_component!,
    protonate_acids, protonate_acids!, deprotonate_oniums, deprotonate_oniums!,
    depolarize, depolarize!, polarize, polarize!,
    to_triple_bond, to_triple_bond!, to_allene_like, to_allene_like!


# TODO: large conjugated system
# TODO: salts and waters should detected by functional group analysis
# TODO: Phosphate, diphosphate, sulfate, nitrate, acetate,
# maleate, fumarate, succinate, citrate, tartrate, oxalate,
# mesylate, tosylate, besylate,
# benzoate, gluconate


"""
    kekulize(mol::SimpleMolGraph) -> Vector{Int}

Kekulize the molecule that has SMILES aromatic bonds.

SMILES allows aromatic atoms in small letters - b, c, n, o, p, s, [as] and [se]. Once these are stored in `SmilesAtom.isaromatic` field, then `kekulize` will place double bonds to satisfy valences.

Kekulization is necessary for molecules parsed from SMILES. If not kekulized, some bond valence and implicit hydrogen properties would be wrong.
"""
function kekulize(mol::SimpleMolGraph{T,V,E}) where {T,V,E}
    nodes = T[]
    for i in vertices(mol)
        get_prop(mol, i, :isaromatic) === nothing && continue
        get_prop(mol, i, :isaromatic) || continue
        if get_prop(mol, i, :symbol) === :C
            push!(nodes, i)
        elseif get_prop(mol, i, :symbol) in (:N, :P, :As)
            if (degree(mol.graph, i) == 2 ||
                    (degree(mol.graph, i) == 3 && get_prop(mol, i, :charge) == 1))
                push!(nodes, i)
            end
        end
    end
    subg, vmap = induced_subgraph(mol.graph, nodes)
    matching = max_matching(subg)
    is_perfect_matching(subg, matching) || error(
        "Kekulization failed: Please check if your SMILES is valid (e.g. Pyrrole n should be [nH])")
    arr = getproperty.(eprops(mol), :order)
    for e in matching
        ge = undirectededge(T, vmap[src(e)], vmap[dst(e)])
        arr[edge_rank(mol, ge)] = 2
    end
    return arr
end
kekulize!(mol::MolGraph) = set_prop!(mol, :e_order, kekulize(mol))


"""
    removable_hydrogens(mol::SimpleMolGraph) -> Set{Int}

Return a vector of removable hydrogen nodes.

removable means not charged, no unpaired electron, no specific mass,
non-stereospecific and attached to organic heavy atoms.
"""
function removable_hydrogens(mol::SimpleMolGraph{T,V,E}) where {T,V,E}
    hs = T[]
    organic_heavy = Set([
        :B, :C, :N, :O, :F, :Si, :P, :S, :Cl, :As, :Se, :Br, :I
    ])
    for i in vertices(mol)
        get_prop(mol, i, :symbol) === :H || continue
        get_prop(mol, i, :charge) == 0 || continue
        get_prop(mol, i, :multiplicity) == 1 || continue
        get_prop(mol, i, :mass) === nothing || continue
        degree(mol.graph, i) == 1 || continue
        nbr = neighbors(mol, i)[1]
        get_prop(mol, nbr, :symbol) in organic_heavy || continue
        get_prop(mol, i, nbr, :order) == 1 || continue
        # TODO: check stereo
        push!(hs, i)
    end
    return hs
end


"""
    allhydrogens(mol::SimpleMolGraph) -> Set{Int}

Return a set of all hydrogen nodes.
"""
function all_hydrogens(mol::SimpleMolGraph{T,V,E}) where {T,V,E}
    hs = T[]
    for i in vertices(mol)
        get_prop(mol, i, :symbol) === :H || continue
        push!(hs, i)
    end
    return hs
end


"""
    removehydrogens(mol::MolGraph) -> GraphMol

Return the molecule with hydrogen nodes removed.

If option `all` is set to true (default), all hydrogens will be removed, otherwise only trivial hydrogens will be removed (see [`trivialhydrogens`](@ref)).
"""
function remove_hydrogens!(mol::SimpleMolGraph; all=true)
    hydrogens = all ? all_hydrogens : removable_hydrogens
    rem_vertices!(mol, hydrogens(mol))
end


"""
    addhydrogens(mol::GraphMol) -> GraphMol

Return the molecule with all hydrogen nodes explicitly attached.
"""
function add_hydrogens!(mol::MolGraph{T,V,E}) where {T,V,E}
    implicit_hs = implicit_hydrogens(mol)
    for i in vertices(mol)
        for j in 1:implicit_hs[i]
            add_vertex!(mol, V(:H))
            add_edge!(mol, i, nv(mol), E())
        end
    end
end



"""
    largest_component_nodes(mol::GraphMol) -> Set{Int}

Return a set of nodes in the largest connected component.
"""
largest_component_nodes(mol::SimpleMolGraph{T,V,E}
    ) where {T,V,E} = sortstablemax(connected_components(mol.graph), by=length, init=T[])


"""
    extractlargestcomponent(mol::GraphMol) -> GraphMol

Return the largest connected component of the molecular graph.

This should be useful when you want to remove salt and water molecules from the molecular graph simply. On the other hand, this can remove important components from the mixture so carefully apply this preprocess method.
"""
extract_largest_component!(mol::MolGraph
    ) = rem_vertices!(mol, setdiff(vertices(mol), largest_component_nodes(mol)))



"""
    protonate_acids(mol::SimpleMolGraph) -> GraphMol

Protonate oxo(thio) anion groups of the molecule.
"""
function protonate_acids(mol::SimpleMolGraph)
    atomsymbol_ = atomsymbol(mol)
    charge_ = charge(mol)
    connectivity_ = connectivity(mol)
    arr = copy(charge_)
    for o in findall(charge_ .== -1)
        atomsymbol_[o] in (:O, :S) || continue
        @assert connectivity_[o] == 1
        nbr = neighbors(mol, o)[1]
        charge_[nbr] == 1 && continue  # polarized double bond
        arr[o] = 0
    end
    return arr
end
protonate_acids!(mol::MolGraph) = set_prop!(mol, :v_charge, protonate_acids(mol))


"""
    deprotonate_oniums(mol::SimpleMolGraph) -> Nothing

Deprotonate onium groups of the molecule.
"""
function deprotonate_oniums(mol::SimpleMolGraph)
    hydrogens_ = total_hydrogens(mol)
    charge_ = charge(mol)
    arr = copy(charge_)
    for o in findall(charge_ .== 1)
        hydrogens_[o] > 0 || continue
        arr[o] = 0
    end
    return arr
end
deprotonate_oniums!(mol::MolGraph) = set_prop!(mol, :v_charge, deprotonate_oniums(mol))


"""
    depolarize(mol::SimpleMolGraph; negative=:O, positive=[:C, :P]) -> Nothing

Depolarize dipole double bonds of the molecule.
"""
function depolarize(mol::SimpleMolGraph; negative=:O, positive=[:C, :P])
    atomsymbol_ = atomsymbol(mol)
    charge_ = charge(mol)
    connectivity_ = connectivity(mol)
    isaromatic_ = isaromatic(mol)
    carr = copy(charge_)
    oarr = copy(bond_order(mol))
    for o in findall(charge_ .== -1)
        atomsymbol_[o] === negative || continue
        @assert connectivity_[o] == 1
        nbr = neighbors(mol, o)[1]
        atomsymbol_[nbr] in positive || continue
        charge_[nbr] == 1 || continue
        isaromatic_[nbr] && continue
        carr[o] = 0
        carr[nbr] = 0
        oarr[edge_rank(mol, o, nbr)] = 2
    end
    return carr, oarr
end
function depolarize!(mol::MolGraph)
    carr, oarr = depolarize(mol)
    set_prop!(mol, :v_charge, carr)
    set_prop!(mol, :e_order, oarr)
    return carr, oarr
end


"""
    polarize(mol::SimpleMolGraph; negative=:O, positive=[:N, :S]) -> Nothing

Polarize dipole double bonds of the molecule.
"""
function polarize(mol::SimpleMolGraph; negative=:O, positive=[:N, :S])
    atomsymbol_ = atomsymbol(mol)
    charge_ = charge(mol)
    bondorder_ = bondorder(mol)
    connectivity_ = connectivity(mol)
    carr = copy(charge_)
    oarr = copy(bondorder_)
    for o in findall(connectivity_ .== 1)
        atomsymbol_[o] === negative || continue
        charge_[o] == 0 || continue
        nbr = neighbors(mol, o)[1]
        atomsymbol_[nbr] in positive || continue
        charge_[nbr] == 0 || continue
        bondorder_[edge_rank(mol, o, nbr)] == 2 || continue
        carr[o] = -1
        carr[nbr] = 1
        oarr[edge_rank(mol, o, nbr)] = 1
    end
    return carr, oarr
end
function polarize!(mol::MolGraph)
    carr, oarr = polarize(mol)
    set_prop!(mol, :v_charge, carr)
    set_prop!(mol, :e_order, oarr)
    return carr, oarr
end



function find_dipoles(mol::SimpleMolGraph)
    charge_ = charge(mol)
    pie_ = apparent_valence(mol) - degree(mol)
    triads = Tuple{Int,Int,Int}[]
    for c in findall((pie_ .== 2) .* (charge_ .== 1))
        negs = Int[]
        notnegs = Int[]
        for nbr in neighbors(mol, c)
            push!(charge_[nbr] == -1 ? negs : notnegs, nbr)
        end
        (length(negs) == 1 && length(notnegs) == 1) || continue
        push!(triads, (negs[1], c, notnegs[1]))
    end
    return triads
end


"""
    to_triple_bond(mol::SimpleMolGraph) -> Nothing

Standardize the molecule so that all 1,3-dipole groups are represented as triple bond and single bond (e.g. Diazo group C=[N+]=[N-] -> [C-][N+]#N).
"""
function to_triple_bond(mol::SimpleMolGraph)
    carr = copy(charge(mol))
    oarr = copy(bond_order(mol))
    for (first, center, third) in find_dipoles(mol)
        get_prop(mol, first, center, :order) == 2 || continue
        carr[first] = 0
        carr[third] = -1
        oarr[edge_rank(mol, first, center)] = 3
        oarr[edge_rank(mol, center, third)] = 1
    end
    return carr, oarr
end

function to_triple_bond!(mol::MolGraph)
    carr, oarr = to_triple_bond(mol)
    set_prop!(mol, :v_charge, carr)
    set_prop!(mol, :e_order, oarr)
    return carr, oarr
end


"""
    to_allene_like(mol::SimpleMolGraph) -> Nothing

Standardize the molecule so that all 1,3-dipole groups are represented as allene-like structure (e.g. Diazo group [C-][N+]#N -> C=[N+]=[N-]).
"""
function to_allene_like(mol::SimpleMolGraph)
    carr = copy(charge(mol))
    oarr = copy(bond_order(mol))
    for (first, center, third) in find_dipoles(mol)
        get_prop(mol, first, center, :order) == 1 || continue
        carr[first] = 0
        carr[third] = -1
        oarr[edge_rank(mol, first, center)] = 2
        oarr[edge_rank(mol, center, third)] = 2
    end
    return carr, oarr
end

function to_allene_like!(mol::MolGraph)
    carr, oarr = to_allene_like(mol)
    set_prop!(mol, :v_charge, carr)
    set_prop!(mol, :e_order, oarr)
    return carr, oarr
end
