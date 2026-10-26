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
    kekulize(mol::SimpleMolGraph{T}) -> Tuple{Vector{Int},Vector{T}}

Return a tuple of the bond orders with kekulization applied and the pyrrole-like aromatic
atoms (e.g. `[nH]`).

Double bonds and single bonds will be assigned to aromatic rings which consist of SMILES
lowercase atoms (called Kekulization). Kekulization is necessary for the valence and
implicit hydrogens of a molecule parsed from SMILES to be correctly evaluated.
An aromatic atom without double bonds gets a double bond if the total bond order (aromatic
bonds as single bonds) is lower than the typical valence (e.g. `c`, `n`, `[nH+]` and `[o+]`,
but not `[nH]`, `o`, `[se]` and `[cH-]`). An aromatic atom with a double bond (e.g. `c=O`
and `[n+](=c)`) gets no more double bonds, except hypervalent atoms (e.g. `n(=O)` of
pyridine N-oxide). A pyridine-like `n` in an odd-membered aromatic system is a pyrrole-like
atom with an implicit hydrogen (e.g. `n1cccc1`).
"""
function kekulize(mol::SimpleMolGraph{T}) where T
    # "raw" bond orders
    bondorder = [bond_order(mol[e]) for e in edges(mol)]
    ernk = edge_rank(mol)
    arom = atom_isaromatic(mol)
    sym_arr = atom_symbol(mol)
    chg_arr = Int[atom_charge(mol[i]) for i in vertices(mol)]
    hs_arr = Int[has_hydrogens(typeof(mol[i])) ? mol[i].hydrogens : 0 for i in vertices(mol)]
    bondsum = apparent_valence(mol.graph, bondorder)
    hasdouble = falses(nv(mol))
    for (e, o) in zip(edges(mol), bondorder)
        o == 2 && (hasdouble[src(e)] = hasdouble[dst(e)] = true)
    end
    # lone pair in p-orbital, pyrrole-like aromatic atom
    pyrrole_like = T[]
    for ring in fused_rings(mol.graph)
        arom_vs = T[]
        canbepyl = T[]  # can be pyrrole-like aromatic atom
        for i in ring
            arom[i] || continue
            if i in mol.gprops.pyrrole_like
                push!(pyrrole_like, i) # explicit pyrrole-like
                continue
            end
            sym = sym_arr[i]
            deg = degree(mol.graph, i)
            deg in (2, 3) || error("Kekulization failed: invalid aromatic $(sym) with valence $(deg)")
            # no pi bond: double bond (except hypervalent) or enough bonds without double bonds
            fv = formal_valence(sym, chg_arr[i])
            if hasdouble[i] ? bondsum[i] <= fv : bondsum[i] + hs_arr[i] >= fv
                # explicit pyrrole-like [nH] (three single bonds including the hydrogen)
                if (sym in (:N, :P, :As) && chg_arr[i] == 0 && deg + hs_arr[i] == 3 && !hasdouble[i]
                        && (hs_arr[i] == 1 || any(sym_arr[nbr] === :H for nbr in neighbors(mol, i))))
                    push!(pyrrole_like, i)
                end
                continue  # o, s, [se], [cH-], c=O, [nH0X3], [n+](=c)
            end
            if deg == 2 && sym in (:N, :P, :As)
                # expected pyridyl, but wrongly can be implicit pyrrole-like
                push!(canbepyl, i)
            end
            push!(arom_vs, i)  # c, [cH0], b, n, n=O, [n+][O-], [o+]
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

"""
    kekulize!(mol::SimpleMolGraph)

In-place version of [`kekulize`](@ref) that stores the bond orders in the descriptors of `mol`
and the pyrrole-like atoms in `mol.gprops.pyrrole_like`. [`apparent_valence`](@ref) is also
overwritten with the total order of the kekulized bonds.
"""
function kekulize!(mol::SimpleMolGraph)
    bondorder, pyrrole_like = kekulize(mol)
    set_descriptor!(mol, Val(:bond_order), bondorder)
    set_descriptor!(mol, Val(:apparent_valence), apparent_valence(mol.graph, bondorder))
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
    all_hydrogens(mol::SimpleMolGraph{T}) -> Vector{T}

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
    remove_hydrogens!(mol::SimpleMolGraph{T}) -> Vector{T}

Remove trivial hydrogen vertices using `Graphs.rem_vertices!` and return the vmap array.

'Trivial' hydrogens are those returned by [`removable_hydrogens`](@ref).
"""
remove_hydrogens!(mol::SimpleMolGraph) = rem_vertices!(mol, removable_hydrogens(mol))


"""
    remove_all_hydrogens!(mol::SimpleMolGraph{T}) -> Vector{T}

Remove all hydrogen vertices using `Graphs.rem_vertices!` and return the vmap array.
"""
remove_all_hydrogens!(mol::SimpleMolGraph) = rem_vertices!(mol, all_hydrogens(mol))


"""
    add_hydrogens!(mol::SimpleMolGraph) -> Nothing

Attach all implicit hydrogens to the molecule as explicit hydrogen vertices.
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
    largest_component_nodes(mol::SimpleMolGraph{T}) -> Vector{T}

Return a vector of nodes in the largest connected component.
"""
largest_component_nodes(mol::SimpleMolGraph{T}
    ) where T = sortstablemax(connected_components(mol.graph), by=length, init=T[])


"""
    extract_largest_component!(mol::SimpleMolGraph{T}) -> Vector{T}

Remove all vertices except for the largest connected component of the molecular graph,
and return the vmap array.

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

"""
    protonate_acids!(mol::SimpleMolGraph)

In-place version of [`protonate_acids`](@ref) that stores the atom charges in the descriptors of `mol`.
"""
function protonate_acids!(mol::SimpleMolGraph)
    set_descriptor!(
        mol, Val(:atom_charge),
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

"""
    deprotonate_oniums!(mol::SimpleMolGraph)

In-place version of [`deprotonate_oniums`](@ref) that stores the atom charges in the descriptors of `mol`.
"""
function deprotonate_oniums!(mol::SimpleMolGraph)
    set_descriptor!(
        mol, Val(:atom_charge),
        deprotonate_oniums(atom_charge(mol), total_hydrogens(mol))
    )
end


"""
    depolarize(mol::SimpleMolGraph; negative=[:O], positive=[:C, :P]
        ) -> Tuple{Vector{Int},Vector{Int}}

Depolarize dipole double bonds of the molecule (e.g. `[C+]-[O-]` -> `C=O`), and return
a tuple of the atom charges and the bond orders.
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

"""
    depolarize!(mol::SimpleMolGraph; negative=[:O], positive=[:C, :P])

In-place version of [`depolarize`](@ref) that stores the atom charges and bond orders
in the descriptors of `mol`.
"""
function depolarize!(mol::SimpleMolGraph; kwargs...)
    carr, oarr = depolarize(
        mol.graph, atom_symbol(mol), atom_charge(mol),
        is_aromatic(mol), bond_order(mol); kwargs...)
    set_descriptor!(mol, Val(:atom_charge), carr)
    set_descriptor!(mol, Val(:bond_order), oarr)
end


"""
    polarize(mol::SimpleMolGraph; negative=[:O], positive=[:N, :S]
        ) -> Tuple{Vector{Int},Vector{Int}}

Polarize dipole double bonds of the molecule (e.g. `N=O` -> `[N+]-[O-]`), and return
a tuple of the atom charges and the bond orders.
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

"""
    polarize!(mol::SimpleMolGraph; negative=[:O], positive=[:N, :S])

In-place version of [`polarize`](@ref) that stores the atom charges and bond orders
in the descriptors of `mol`.
"""
function polarize!(mol::SimpleMolGraph; kwargs...)
    carr, oarr = polarize(
        mol.graph, atom_symbol(mol), atom_charge(mol), bond_order(mol); kwargs...)
    set_descriptor!(mol, Val(:atom_charge), carr)
    set_descriptor!(mol, Val(:bond_order), oarr)
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
    to_triple_bond(mol::SimpleMolGraph) -> Tuple{Vector{Int},Vector{Int}}

Return a tuple of the atom charges and the bond orders where all 1,3-dipole groups are
represented as triple bond and single bond (e.g. Diazo group C=[N+]=[N-] -> [C-][N+]#N).
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

"""
    to_triple_bond!(mol::SimpleMolGraph)

In-place version of [`to_triple_bond`](@ref) that stores the atom charges and bond orders
in the descriptors of `mol`.
"""
function to_triple_bond!(mol::SimpleMolGraph)
    carr, oarr = to_triple_bond(
        mol.graph, atom_charge(mol), apparent_valence(mol), bond_order(mol))
    set_descriptor!(mol, Val(:atom_charge), carr)
    set_descriptor!(mol, Val(:bond_order), oarr)
end


"""
    to_allene_like(mol::SimpleMolGraph) -> Tuple{Vector{Int},Vector{Int}}

Return a tuple of the atom charges and the bond orders where all 1,3-dipole groups are
represented as allene-like structure (e.g. Diazo group [C-][N+]#N -> C=[N+]=[N-]).
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

"""
    to_allene_like!(mol::SimpleMolGraph)

In-place version of [`to_allene_like`](@ref) that stores the atom charges and bond orders
in the descriptors of `mol`.
"""
function to_allene_like!(mol::SimpleMolGraph)
    carr, oarr = to_allene_like(
        mol.graph, atom_charge(mol), apparent_valence(mol), bond_order(mol))
    set_descriptor!(mol, Val(:atom_charge), carr)
    set_descriptor!(mol, Val(:bond_order), oarr)
end
