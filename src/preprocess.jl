#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

function remap!(::Val{:pyrrole_like}, gprop::SimpleMolProperty{T},
        vmap::Vector{T}, edges::Vector{Edge{T}}) where T <: Integer
    revv = Dict(v => i for (i, v) in enumerate(vmap))
    vec = T[revv[v] for v in gprop.pyrrole_like if haskey(revv, v)]
    empty!(gprop.pyrrole_like)
    append!(gprop.pyrrole_like, sort(vec))
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
    bondorder = [bond_order(mol[e]) for e in edges(mol)]
    ernk = edge_rank(mol)
    # lone pair in p-orbital, pyrrole-like aromatic atom
    pyrrole_like = T[]
    for ring in fused_rings(mol.graph)
        arom_vs = T[]
        canbepyl = T[]  # can be pyrrole-like aromatic atom
        for i in ring
            atmp = mol[i]
            has_isaromatic(typeof(atmp)) && isaromatic(atmp) || continue  # not nothing or false
            if i in mol.gprops.pyrrole_like
                push!(pyrrole_like, i) # explicit pyrrole-like
                continue
            end
            sym = atom_symbol(atmp)
            deg = degree(mol.graph, i)
            if deg == 2
                sym in (:O, :S) && continue  # o, s, se
                if sym in (:N, :P, :As)
                    # expected pyridyl, but wrongly can be implicit pyrrole-like
                    push!(canbepyl, i)
                end
                push!(arom_vs, i)  # including [cH0], b
            elseif deg == 3
                hasdouble = any(bondorder[edge_rank(ernk, i, nbr)] == 2 for nbr in neighbors(mol, i))
                if sym === :C
                    hasdouble && continue  # c=O
                    push!(arom_vs, i)
                else  # sym in (:N, :P, :As)
                    if atom_charge(atmp) == 0 && !hasdouble
                        # explicit pyrrole-like [nH]
                        if has_hydrogens(typeof(atmp)) && atmp.hydrogens == 1
                            push!(pyrrole_like, i)
                        elseif any(atom_symbol(mol[nbr]) === :H for nbr in neighbors(mol, i))
                            push!(pyrrole_like, i)
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
            bondorder[edge_rank(ernk, vmap[src(e)], vmap[dst(e)])] = 2
        end
    end
    return bondorder, pyrrole_like
end

function kekulize!(mol::SimpleMolGraph)
    bondorder, pyrrole_like = kekulize(mol)
    set_descriptor!(mol, :bond_order, bondorder)
    empty!(mol.gprops.pyrrole_like)
    append!(mol.gprops.pyrrole_like, pyrrole_like)
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
        atom_symbol(mol[i]) === :H || continue
        atom_charge(mol[i]) == 0 || continue
        multiplicity(mol[i]) == 1 || continue
        isotope(mol[i]) == 0 || continue
        degree(mol, i) == 1 || continue
        nbr = only(neighbors(mol, i))
        atom_symbol(mol[nbr]) in organic_heavy || continue
        bond_order(mol[i, nbr]) == 1 || continue
        haskey(mol[:stereocenter], nbr) && continue
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
        atom_symbol(mol[i]) === :H || continue
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
remove_all_hydrogens!(mol::SimpleMolGraph) = rem_vertices!(mol, all_hydrogens(mol))


"""
    add_hydrogens!(mol::SimpleMolGraph) -> Nothing

Return the molecule with all hydrogen nodes explicitly attached.
"""
add_hydrogens!(mol::SimpleMolGraph) = add_hydrogens!(mol, vproptype(mol), eproptype(mol))

function add_hydrogens!(
        mol::SimpleMolGraph, V::Type{<:StandardAtom}, E::Type{<:StandardBond})
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
protonate_acids(mol::SimpleMolGraph
    ) = protonate_acids(mol.graph, atom_symbol(mol), atom_charge(mol))

function protonate_acids(
        g::SimpleGraph, symbol_arr::Vector{Symbol}, charge_arr::Vector{Int})
    new_charge = copy(charge_arr)
    for o in findall(charge_arr .== -1)
        symbol_arr[o] in (:O, :S) || continue
        degree(g, o) == 1 || continue  # atypical valence
        charge_arr[only(neighbors(g, o))] == 1 && continue  # polarized double bond
        new_charge[o] = 0
    end
    return new_charge
end

function protonate_acids!(mol::SimpleMolGraph)
    set_descriptor!(
        mol, :atom_charge,
        protonate_acids(mol.graph, atom_symbol(mol), atom_charge(mol))
    )
end


"""
    deprotonate_oniums(mol::SimpleMolGraph) -> Vector{Int}

Deprotonate onium groups of the molecule.
"""
deprotonate_oniums(mol::SimpleMolGraph
    ) = deprotonate_oniums(atom_charge(mol), total_hydrogens(mol))

function deprotonate_oniums(charge_arr::Vector{Int}, h_arr::Vector{Int})
    new_charge = copy(charge_arr)
    for o in findall(charge_arr .== 1)
        h_arr[o] > 0 || continue
        new_charge[o] = 0
    end
    return new_charge
end

function deprotonate_oniums!(mol::SimpleMolGraph)
    set_descriptor!(
        mol, :atom_charge,
        deprotonate_oniums(atom_charge(mol), total_hydrogens(mol))
    )
end


"""
    depolarize(mol::SimpleMolGraph; negative=:O, positive=[:C, :P]) -> Nothing

Depolarize dipole double bonds of the molecule.
"""
depolarize(mol::SimpleMolGraph; kwargs...) = depolarize(
    mol.graph, atom_symbol(mol), atom_charge(mol),
    is_aromatic(mol), bond_order(mol); kwargs...
)

function depolarize(
        g::SimpleGraph, symbol_arr::Vector{Symbol}, charge_arr::Vector{Int},
        is_arom_arr::BitVector, order_arr::Vector{Int}; negative=[:O], positive=[:C, :P])
    carr = copy(charge_arr)
    oarr = copy(order_arr)
    ernk = edge_rank(g)
    for o in findall(charge_arr .== -1)
        symbol_arr[o] in negative || continue
        degree(g, o) == 1 || continue  # atypical valence
        nbr = only(neighbors(g, o))
        symbol_arr[nbr] in positive || continue
        charge_arr[nbr] == 1 || continue
        is_arom_arr[nbr] && continue
        carr[o] = 0
        carr[nbr] = 0
        oarr[edge_rank(ernk, o, nbr)] = 2
    end
    return carr, oarr
end

function depolarize!(mol::SimpleMolGraph; kwargs...)
    carr, oarr = depolarize(
        mol.graph, atom_symbol(mol), atom_charge(mol),
        is_aromatic(mol), bond_order(mol); kwargs...)
    set_descriptor!(mol, :atom_charge, carr)
    set_descriptor!(mol, :bond_order, oarr)
end


"""
    polarize(mol::SimpleMolGraph; negative=:O, positive=[:N, :S]) -> Nothing

Polarize dipole double bonds of the molecule.
"""
polarize(mol::SimpleMolGraph; kwargs...) = polarize(
    mol.graph, atom_symbol(mol), atom_charge(mol), bond_order(mol); kwargs...
)
function polarize(
        g::SimpleGraph, symbol_arr::Vector{Symbol}, charge_arr::Vector{Int},
        order_arr::Vector{Int}; negative=[:O], positive=[:N, :S])
    carr = copy(charge_arr)
    oarr = copy(order_arr)
    ernk = edge_rank(g)
    for o in findall(degree(g) .== 1)
        symbol_arr[o] in negative || continue
        charge_arr[o] == 0 || continue
        nbr = only(neighbors(g, o))
        symbol_arr[nbr] in positive || continue
        charge_arr[nbr] == 0 || continue
        order_arr[edge_rank(ernk, o, nbr)] == 2 || continue
        carr[o] = -1
        carr[nbr] = 1
        oarr[edge_rank(ernk, o, nbr)] = 1
    end
    return carr, oarr
end

function polarize!(mol::SimpleMolGraph; kwargs...)
    carr, oarr = polarize(
        mol.graph, atom_symbol(mol), atom_charge(mol), bond_order(mol); kwargs...)
    set_descriptor!(mol, :atom_charge, carr)
    set_descriptor!(mol, :bond_order, oarr)
end


function find_dipoles(g::SimpleGraph, charge_arr::Vector{Int}, appval_arr::Vector{Int})
    pie_arr = appval_arr - degree(g)
    triads = Tuple{Int,Int,Int}[]
    for c in findall((pie_arr .== 2) .* (charge_arr .== 1))
        negs = Int[]
        notnegs = Int[]
        for nbr in neighbors(g, c)
            push!(charge_arr[nbr] == -1 ? negs : notnegs, nbr)
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
to_triple_bond(mol::SimpleMolGraph) = to_triple_bond(
    mol.graph, atom_charge(mol), apparent_valence(mol), bond_order(mol))

function to_triple_bond(
        g::SimpleGraph, charge_arr::Vector{Int}, appval_arr::Vector{Int}, order_arr::Vector{Int})
    carr = copy(charge_arr)
    oarr = copy(order_arr)
    ernk = edge_rank(g)
    for (first, center, third) in find_dipoles(g, charge_arr, appval_arr)
        order_arr[edge_rank(ernk, first, center)] == 2 || continue
        carr[first] = 0
        carr[third] = -1
        oarr[edge_rank(ernk, first, center)] = 3
        oarr[edge_rank(ernk, center, third)] = 1
    end
    return carr, oarr
end

function to_triple_bond!(mol::SimpleMolGraph)
    carr, oarr = to_triple_bond(
        mol.graph, atom_charge(mol), apparent_valence(mol), bond_order(mol))
    set_descriptor!(mol, :atom_charge, carr)
    set_descriptor!(mol, :bond_order, oarr)
end


"""
    to_allene_like(mol::SimpleMolGraph) -> Nothing

Standardize the molecule so that all 1,3-dipole groups are represented as allene-like structure (e.g. Diazo group [C-][N+]#N -> C=[N+]=[N-]).
"""
to_allene_like(mol::SimpleMolGraph) = to_allene_like(
    mol.graph, atom_charge(mol), apparent_valence(mol), bond_order(mol))

function to_allene_like(
        g::SimpleGraph, charge_arr::Vector{Int}, appval_arr::Vector{Int}, order_arr::Vector{Int})
    carr = copy(charge_arr)
    oarr = copy(order_arr)
    ernk = edge_rank(g)
    for (first, center, third) in find_dipoles(g, charge_arr, appval_arr)
        order_arr[edge_rank(ernk, first, center)] == 1 || continue
        carr[first] = 0
        carr[third] = -1
        oarr[edge_rank(ernk, first, center)] = 2
        oarr[edge_rank(ernk, center, third)] = 2
    end
    return carr, oarr
end

function to_allene_like!(mol::SimpleMolGraph)
    carr, oarr = to_allene_like(
        mol.graph, atom_charge(mol), apparent_valence(mol), bond_order(mol))
    set_descriptor!(mol, :atom_charge, carr)
    set_descriptor!(mol, :bond_order, oarr)
end
