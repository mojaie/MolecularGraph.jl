#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    kekulize, kekulize!,
    removable_hydrogens, all_hydrogens,
    remove_hydrogens!, remove_all_hydrogens!, add_hydrogens!,
    largest_component_nodes, extract_largest_component!,
    protonate_acids, protonate_acids!, deprotonate_oniums, deprotonate_oniums!,
    depolarize, depolarize!, polarize, polarize!,
    to_triple_bond, to_triple_bond!, to_allene_like, to_allene_like!


struct PyrroleLike{T} <: AbstractVector{T}
    vertices::Vector{T}
end

Base.size(p::PyrroleLike) = size(p.vertices)
Base.getindex(p::PyrroleLike, i...) = getindex(p.vertices, i...)
to_dict(p::PyrroleLike) = p.vertices

remap(p::PyrroleLike{T}, vmap::Dict
    ) where T = PyrroleLike{T}([vmap[v] for v in p.vertices if v in keys(vmap)])


"""
    kekulize(mol::SimpleMolGraph) -> Vector{Int}

Return an array of bond orders with kekulization applied.

Double bonds and single bonds will be assigned to aromatic rings which consist of SMILES
lowercase atoms (called Kekulization). Kekulization is necessary for the valence and
implicit hydrogens of a molecule parsed from SMILES to be correctly evaluated.
"""
function kekulize(mol::SimpleMolGraph{T,V,E}) where {T,V,E}
    bondorder = [get_prop(mol, e, :order) for e in edges(mol)]
    # lone pair in p-orbital, pyrrole-like aromatic atom
    pyrrole_like = T[]
    for ring in fused_rings(mol.graph)
        arom_vs = T[]
        canbepyl = T[]  # can be pyrrole-like aromatic atom
        for i in ring
            get_prop(mol, i, :isaromatic) === true || continue  # not nothing or false
            if has_prop(mol, :pyrrole_like) && i in get_prop(mol, :pyrrole_like)
                push!(pyrrole_like, i)
                continue
            end
            sym = get_prop(mol, i, :symbol)
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
                    if get_prop(mol, i, :charge) == 0 && !hasdouble
                        if any(get_prop(mol, nbr, :symbol) === :H for nbr in neighbors(mol, i))
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
            push!(pyrrole_like, only(canbepyl))
        end
        subg, vmap = induced_subgraph(mol.graph, arom_vs)
        matching = max_matching(subg)
        is_perfect_matching(subg, matching) || error("Kekulization failed: invalid aromatic ring")
        for e in matching
            bondorder[edge_rank(mol, vmap[src(e)], vmap[dst(e)])] = 2
        end
    end
    return bondorder, PyrroleLike(pyrrole_like)
end

function kekulize!(mol::MolGraph)
    bondorder, pyrrole_like = kekulize(mol)
    set_cache!(mol, :e_order, bondorder)
    mol.gprops[:pyrrole_like] = pyrrole_like
end


"""
    removable_hydrogens(mol::SimpleMolGraph{T,V,E}) -> Vector{T}

Return a vector of removable hydrogen nodes.

Removable hydrogens are not charged, have no unpaired electron, have no specific mass,
are non-stereospecific and are attached to organic heavy atoms.
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
        haskey(get_prop(mol, :stereocenter), nbr) && continue
        push!(hs, i)
    end
    return hs
end


"""
    all_hydrogens(mol::SimpleMolGraph{T,V,E}) -> Vector{T}

Return a vector of all hydrogen nodes.
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
    remove_hydrogens!(mol::SimpleMolGraph{T,V,E}) -> Vector{T}

Remove following hydrogen vertices from the molecule: that are not charged, have no
unpaired electron, have no specific mass, are non-stereospecific and are
attached to organic heavy atoms.

This returns vmap array similar to `Graphs.rem_vertices!`.
"""
remove_hydrogens!(mol::SimpleMolGraph) = rem_vertices!(mol, removable_hydrogens(mol))


"""
    remove_all_hydrogens!(mol::SimpleMolGraph{T,V,E}) -> Vector{T}

Remove all hydrogen vertices from the molecule.

This returns vmap array similar to `Graphs.rem_vertices!`.
"""
function remove_all_hydrogens!(mol::SimpleMolGraph{T,V,E}) where {T,V,E}
    to_remove = T[]
    for center in keys(get_prop(mol, :stereocenter))
        safe_stereo_hydrogen!(mol, center)
    end
    vmap = rem_vertices!(mol, all_hydrogens(mol))
    return vmap
end


"""
    add_hydrogens!(mol::SimpleMolGraph)

Return the molecule with all hydrogen nodes explicitly attached.
"""
function add_hydrogens!(mol::SimpleMolGraph{T,V,E}) where {T,V,E}
    implicit_hs = implicit_hydrogens(mol)
    for i in vertices(mol)
        for j in 1:implicit_hs[i]
            add_vertex!(mol, V(:H))
            add_edge!(mol, i, nv(mol), E())
        end
    end
end



"""
    largest_component_nodes(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of nodes in the largest connected component.
"""
largest_component_nodes(mol::SimpleMolGraph{T,V,E}
    ) where {T,V,E} = sortstablemax(connected_components(mol.graph), by=length, init=T[])


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

protonate_acids!(mol::MolGraph) = set_cache!(mol, :v_charge, protonate_acids(mol))


"""
    deprotonate_oniums(mol::SimpleMolGraph) -> Vector{Int}

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

deprotonate_oniums!(mol::MolGraph) = set_cache!(mol, :v_charge, deprotonate_oniums(mol))


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
    set_cache!(mol, :v_charge, carr)
    set_cache!(mol, :e_order, oarr)
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
    set_cache!(mol, :v_charge, carr)
    set_cache!(mol, :e_order, oarr)
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
    set_cache!(mol, :v_charge, carr)
    set_cache!(mol, :e_order, oarr)
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
    set_cache!(mol, :v_charge, carr)
    set_cache!(mol, :e_order, oarr)
    return carr, oarr
end
