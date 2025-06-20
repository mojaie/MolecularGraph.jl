#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

function remap!(::Val{:pyrrole_like}, gprop::SimpleMolProperty{T},
        vmap::Vector{T}, edges::Vector{Edge{T}}) where T <: Integer
    revv = Dict(v => i for (i, v) in enumerate(vmap))
    vec = T[revv[v] for v in gprop.pyrrole_like if haskey(revv, v)]
    empty!(gprop.pyrrole_like)
    append!(gprop.pyrrole_like, vec)
    return
end


"""
    kekulize(mol::SimpleMolGraph) -> Vector{Int}

Return an array of bond orders with kekulization applied.

Double bonds and single bonds will be assigned to aromatic rings which consist of SMILES
lowercase atoms (called Kekulization). Kekulization is necessary for the valence and
implicit hydrogens of a molecule parsed from SMILES to be correctly evaluated.
"""
function kekulize(mol::SimpleMolGraph{T}) where T
    # "raw" bond orders
    bondorder = [bond_order(props(mol, e)) for e in edges(mol)]
    # lone pair in p-orbital, pyrrole-like aromatic atom
    pyrrole_like = T[]
    for ring in fused_rings(mol.graph)
        arom_vs = T[]
        canbepyl = T[]  # can be pyrrole-like aromatic atom
        for i in ring
            get_prop(mol, i, :isaromatic) === true || continue  # not nothing or false
            if i in mol.gprops.pyrrole_like
                push!(pyrrole_like, i) # explicit pyrrole-like
                continue
            end
            sym = atom_symbol(props(mol, i))
            deg = degree(mol.graph, i)
            if deg == 2
                sym in (:O, :S) && continue  # o, s, se
                if sym in (:N, :P, :As)
                    # expected pyridyl, but wrongly can be implicit pyrrole-like
                    push!(canbepyl, i)
                end
                push!(arom_vs, i)  # including [cH0], b
            elseif deg == 3
                hasdouble = any(bondorder[edge_rank(mol, i, nbr)] == 2 for nbr in neighbors(mol, i))
                if sym === :C
                    hasdouble && continue  # c=O
                    push!(arom_vs, i)
                else  # sym in (:N, :P, :As)
                    if atom_charge(props(mol, i)) == 0 && !hasdouble
                        if any(atom_symbol(props(mol, nbr)) === :H for nbr in neighbors(mol, i))
                            push!(pyrrole_like, i) # explicit pyrrole-like [nH]
                        end
                        continue  # including [nH0X3]
                    end
                    push!(arom_vs, i)  # n=O, [n+][O-]
                end
            else
                error("Kekulization failed: invalid aromatic $(sym) with valence $(deg)")
            end
        end
        if length(arom_vs) % 2 == 1
            isempty(canbepyl) && error("Kekulization failed: invalid aromatic ring")
            length(canbepyl) > 1 && error("Kekulization failed: ambiguity in pyrrole-like ring")
            # if there is only one n without hydrogen, that might be a pyrrole-like n
            setdiff!(arom_vs, canbepyl)
            # push!(pyrrole_like, only(canbepyl))  # implicit pyrrole n is not stored
        end
        subg, vmap = induced_subgraph(mol.graph, arom_vs)
        matching = max_matching(subg)
        is_perfect_matching(subg, matching) || error("Kekulization failed: invalid aromatic ring")
        for e in matching
            bondorder[edge_rank(mol, vmap[src(e)], vmap[dst(e)])] = 2
        end
    end
    return bondorder, pyrrole_like
end

function kekulize!(mol::ReactiveMolGraph)
    bondorder, pyrrole_like = kekulize(mol)
    mol.gprops.descriptors.bond_order = bondorder
    mol.gprops.pyrrole_like = pyrrole_like
    return
end


"""
    removable_hydrogens(mol::SimpleMolGraph{T}) -> Vector{T}

Return a vector of removable hydrogen nodes.

Removable hydrogens are not charged, have no unpaired electron, have no specific mass,
are non-stereospecific and are attached to organic heavy atoms.
"""
function removable_hydrogens(mol::SimpleMolGraph{T}) where T
    hs = T[]
    organic_heavy = Set([
        :B, :C, :N, :O, :F, :Si, :P, :S, :Cl, :As, :Se, :Br, :I
    ])
    for i in vertices(mol)
        atom_symbol(props(mol, i)) === :H || continue
        atom_charge(props(mol, i)) == 0 || continue
        multiplicity(props(mol, i)) == 1 || continue
        isnothing(atom_mass(props(mol, i))) || continue
        degree(mol.graph, i) == 1 || continue
        nbr = only(neighbors(mol, i))
        atom_symbol(props(mol, nbr)) in organic_heavy || continue
        bond_order(props(mol, i, nbr)) == 1 || continue
        haskey(get_prop(mol, :stereocenter), nbr) && continue
        push!(hs, i)
    end
    return hs
end


"""
    all_hydrogens(mol::SimpleMolGraph) -> Vector{T}

Return a vector of all hydrogen nodes.
"""
function all_hydrogens(mol::SimpleMolGraph{T}) where T
    hs = T[]
    for i in vertices(mol)
        atom_symbol(props(mol, i)) === :H || continue
        push!(hs, i)
    end
    return hs
end


"""
    remove_hydrogens!(mol::SimpleMolGraph) -> Vector{T<:Integer}

Remove trivial hydrogen vertices using `Graphs.rem_vertices!` and return vmap array.

'Trivial' means that are not charged, have no unpaired electron, have no specific mass,
are non-stereospecific and are attached to organic heavy atoms.
"""
remove_hydrogens!(mol::SimpleMolGraph) = rem_vertices!(mol, removable_hydrogens(mol))


"""
    remove_all_hydrogens!(mol::SimpleMolGraph) -> Vector{T<:Integer}

Remove all hydrogen vertices using `Graphs.rem_vertices!` and return vmap array.
"""
function remove_all_hydrogens!(mol::SimpleMolGraph)
    for center in keys(mol.gprops.stereocenter)
        safe_stereo_hydrogen!(mol, center)
    end
    vmap = rem_vertices!(mol, all_hydrogens(mol))
    return vmap
end


"""
    add_hydrogens!(mol::ReactiveMolGraph) -> Nothing

Return the molecule with all hydrogen nodes explicitly attached.
"""
function add_hydrogens!(mol::ReactiveMolGraph{T,V,E}) where {T,V,E}
    implicit_hs = implicit_hydrogens(mol)
    for i in vertices(mol)
        for j in 1:implicit_hs[i]
            add_vertex!(mol, V(;symbol=:H))
            add_edge!(mol, i, nv(mol), E())
        end
    end
    return
end



"""
    largest_component_nodes(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of nodes in the largest connected component.
"""
largest_component_nodes(mol::SimpleMolGraph{T}
    ) where T = sortstablemax(connected_components(mol.graph), by=length, init=T[])


"""
    extract_largest_component!(mol::SimpleMolGraph) -> Nothing

Return the largest connected component of the molecular graph.

This should be useful when you want to remove salt and water molecules from the molecular graph simply. On the other hand, this can remove important components from the mixture so carefully apply this preprocess method.
"""
extract_largest_component!(mol::SimpleMolGraph
    ) = rem_vertices!(mol, setdiff(vertices(mol), largest_component_nodes(mol)))



"""
    protonate_acids(mol::SimpleMolGraph) -> Vector{Int}

Protonate oxo(thio) anion groups of the molecule.
"""
function protonate_acids(mol::SimpleMolGraph)
    atomsymbol_ = atom_symbol(mol)
    charge_ = atom_charge(mol)
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

function protonate_acids!(mol::ReactiveMolGraph)
    has_descriptor(mol, :atom_charge) || error("Descriptor :atom_charge not available")
    set_descriptor!(mol, :atom_charge, protonate_acids(mol))
    return
end

"""
    deprotonate_oniums(mol::SimpleMolGraph) -> Vector{Int}

Deprotonate onium groups of the molecule.
"""
function deprotonate_oniums(mol::SimpleMolGraph)
    hydrogens_ = total_hydrogens(mol)
    charge_ = atom_charge(mol)
    arr = copy(charge_)
    for o in findall(charge_ .== 1)
        hydrogens_[o] > 0 || continue
        arr[o] = 0
    end
    return arr
end

function deprotonate_oniums!(mol::ReactiveMolGraph)
    has_descriptor(mol, :atom_charge) || error("Descriptor :atom_charge not available")
    set_descriptor!(mol, :atom_charge, deprotonate_oniums(mol))
    return
end


"""
    depolarize(mol::SimpleMolGraph; negative=:O, positive=[:C, :P]) -> Nothing

Depolarize dipole double bonds of the molecule.
"""
function depolarize(mol::SimpleMolGraph; negative=:O, positive=[:C, :P])
    atomsymbol_ = atom_symbol(mol)
    charge_ = atom_charge(mol)
    connectivity_ = connectivity(mol)
    isaromatic_ = is_aromatic(mol)
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

function depolarize!(mol::ReactiveMolGraph)
    has_descriptor(mol, :atom_charge) || error("Descriptor :atom_charge not available")
    has_descriptor(mol, :bond_order) || error("Descriptor :bond_order not available")
    carr, oarr = depolarize(mol)
    set_descriptor!(mol, :atom_charge, carr)
    set_descriptor!(mol, :bond_order, oarr)
    return
end


"""
    polarize(mol::SimpleMolGraph; negative=:O, positive=[:N, :S]) -> Nothing

Polarize dipole double bonds of the molecule.
"""
function polarize(mol::SimpleMolGraph; negative=:O, positive=[:N, :S])
    atomsymbol_ = atom_symbol(mol)
    charge_ = atom_charge(mol)
    bondorder_ = bond_order(mol)
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

function polarize!(mol::ReactiveMolGraph)
    has_descriptor(mol, :atom_charge) || error("Descriptor :atom_charge not available")
    has_descriptor(mol, :bond_order) || error("Descriptor :bond_order not available")
    carr, oarr = polarize(mol)
    set_descriptor!(mol, :atom_charge, carr)
    set_descriptor!(mol, :bond_order, oarr)
    return
end


function find_dipoles(mol::SimpleMolGraph)
    charge_ = atom_charge(mol)
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
    carr = copy(atom_charge(mol))
    oarr = copy(bond_order(mol))
    for (first, center, third) in find_dipoles(mol)
        bond_order(props(mol, first, center)) == 2 || continue
        carr[first] = 0
        carr[third] = -1
        oarr[edge_rank(mol, first, center)] = 3
        oarr[edge_rank(mol, center, third)] = 1
    end
    return carr, oarr
end

function to_triple_bond!(mol::ReactiveMolGraph)
    has_descriptor(mol, :atom_charge) || error("Descriptor :atom_charge not available")
    has_descriptor(mol, :bond_order) || error("Descriptor :bond_order not available")
    carr, oarr = to_triple_bond(mol)
    set_descriptor!(mol, :atom_charge, carr)
    set_descriptor!(mol, :bond_order, oarr)
    return
end

"""
    to_allene_like(mol::SimpleMolGraph) -> Nothing

Standardize the molecule so that all 1,3-dipole groups are represented as allene-like structure (e.g. Diazo group [C-][N+]#N -> C=[N+]=[N-]).
"""
function to_allene_like(mol::SimpleMolGraph)
    carr = copy(atom_charge(mol))
    oarr = copy(bond_order(mol))
    for (first, center, third) in find_dipoles(mol)
        bond_order(props(mol, first, center)) == 1 || continue
        carr[first] = 0
        carr[third] = -1
        oarr[edge_rank(mol, first, center)] = 2
        oarr[edge_rank(mol, center, third)] = 2
    end
    return carr, oarr
end

function to_allene_like!(mol::ReactiveMolGraph)
    has_descriptor(mol, :atom_charge) || error("Descriptor :atom_charge not available")
    has_descriptor(mol, :bond_order) || error("Descriptor :bond_order not available")
    carr, oarr = to_allene_like(mol)
    set_descriptor!(mol, :atom_charge, carr)
    set_descriptor!(mol, :bond_order, oarr)
    return
end