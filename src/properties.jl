#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    sssr, sssr!,
    which_ring, edge_which_ring, fused_rings, which_fused_ring,
    smallest_ring, ring_count, is_in_ring, is_edge_in_ring,
    atom_symbol, charge, multiplicity, bond_order,
    lone_pair, lone_pair!, apparent_valence, apparent_valence!, valence, valence!,
    explicit_hydrogens, implicit_hydrogens, heavy_atoms,
    total_hydrogens, connectivity,
    is_hydrogen_donor, hydrogen_donor_count,
    is_hydrogen_acceptor, hydrogen_acceptor_count,
    is_rotatable, rotatable_count,
    atom_counter, heavy_atom_count, molecular_formula, empirical_formula,
    pi_electron, pi_delocalized, hybridization, hybridization_delocalized,
    is_ring_aromatic, is_ring_aromatic!, is_aromatic, is_edge_aromatic


const LONEPAIR_COUNT = Dict(
    :H => 0, :B => -1, :C => 0, :N => 1, :O => 2, :F => 3,
    :Si => 0, :P => 1, :S => 2, :Cl => 3,
    :As => 1, :Se => 2, :Br => 3, :I => 3
)



# Molecular graph topology descriptors

"""
    degree(mol::SimpleMolGraph{T}) -> Vector{T}

Return a vector of size ``n`` representing the node degree of the molecular graph
of 1 to ``n``th atoms of the given molecule.

This property corresponds to SMARTS `D` query.
"""
function Graphs.degree(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_degree) && return get_cache(mol, :v_degree)
    return degree(mol.graph)
end

degree!(mol::SimpleMolGraph) = set_cache!(mol, :v_degree, degree(mol.graph))


"""
    sssr(mol::SimpleMolGraph{T}) -> Vector{Vector{T}}

Return vectors of ring nodes representing small set of smallest rings (SSSR).

See [`mincyclebasis`](@ref).
"""
function sssr(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :sssr) && return get_cache(mol, :sssr)
    return mincyclebasis(mol.graph)
end

sssr!(mol::SimpleMolGraph) = set_cache!(mol, :sssr, mincyclebasis(mol.graph))


"""
    which_ring(mol::SimpleMolGraph) -> Vector{Vector{Int}}

Return a vector of size ``n`` representing [`sssr`](@ref) membership of
1 to ``n``th nodes of the given graph.

SSSR membership is represented as a vector of SSSR indices assigned to each rings.
This means nodes that have the same SSSR index belong to the same SSSR.
"""
function which_ring(mol::SimpleMolGraph{T,V,E}) where {T,V,E}
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_which_ring) && return get_cache(mol, :v_which_ring)
    arr = [Int[] for _ in vertices(mol)]
    for (i, cyc) in enumerate(sssr(mol))
        for n in cyc
            push!(arr[n], i)
        end
    end
    return arr
end


"""
    edge_which_sssr(mol::SimpleMolGraph) -> Vector{Vector{Int}}

Return a vector of size ``n`` representing [`sssr`](@ref) membership of
1 to ``n``th bonds of the given molecule.

SSSR membership is represented as a set of SSSR indices assigned to each rings.
This means bonds that have the same SSSR index belong to the same SSSR.
"""
function edge_which_ring(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :e_which_ring) && return get_cache(mol, :e_which_ring)
    arr = [Int[] for _ in 1:ne(mol)]
    for (i, cyc) in enumerate(sssr(mol))
        for j in 1:(length(cyc) - 1)
            push!(arr[edge_rank(mol, cyc[j], cyc[j + 1])], i)
        end
        push!(arr[edge_rank(mol, cyc[1], cyc[end])], i)
    end
    return arr
end


"""
    fused_rings(mol::SimpleMolGraph{T}) -> Vector{Vector{T}}

Return vectors of fused ring node sets.

A fused ring is defined as a 2-edge connected components in terms of graph theory.
Spirocyclic structures are considered to be part of a fused ring.
"""
function fused_rings(g::SimpleGraph)
    cobr = setdiff(Set(edges(g)), bridges(g))
    subg, vmap = induced_subgraph(g, collect(cobr))
    return  [vmap[c] for c in connected_components(subg)]
end

function fused_rings(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :fused_rings) && return get_cache(mol, :fused_rings)
    return fused_rings(mol.graph)
end


"""
    which_fused_ring(mol::SimpleMolGraph) -> Vector{Vector{Int}}

Return a vector of size ``n`` representing [`fused_rings`](@ref) membership of
1 to ``n``th atoms of the given molecule.

Fused ring membership is represented as a set of fused ring indices assigned to each fused rings.
This means atoms that have the same fused ring index belong to the same fused ring.
"""
function which_fused_ring(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_which_fused_ring) && return get_cache(mol, :v_which_fused_ring)
    arr = [Int[] for _ in vertices(mol)]
    for (i, conn) in enumerate(fused_rings(mol))
        for n in conn
            push!(arr[n], i)
        end
    end
    return arr
end


"""
    smallest_ring(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the size of the smallest [`sssr`](@ref)
that 1 to ``n``th atoms of the given molecule belong to. 

If the node is not in a ring, the value would be 0.
This property corresponds to SMARTS `r` query.
"""
function smallest_ring(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_smallest_ring) && return get_cache(mol, :v_smallest_ring)
    sssr_ = sssr(mol)
    whichring_ = which_ring(mol)
    arr = zeros(Int, nv(mol))
    for i in vertices(mol)
        rs = whichring_[i]
        isempty(rs) && continue
        arr[i] = minimum(length, sssr_[rs])
    end
    return arr
end


"""
    ring_count(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of [`sssr`](@ref)
that 1 to ``n``th atoms of the given molecule belong to.

This property corresponds to SMARTS `R` query.
"""
function ring_count(mol::MolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_ring_count) && return get_cache(mol, :v_ring_count)
    return length.(which_ring(mol))
end


"""
    is_in_ring(mol::MolGraph) -> Vector{Bool}

Return a vector of size ``n`` representing whether 1 to ``n``th atoms of
the given molecule belong to a ring or not.
"""
function is_in_ring(mol::MolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_is_in_ring) && return get_cache(mol, :v_is_in_ring)
    return .!isempty.(which_ring(mol))
end


"""
    is_edge_in_ring(mol::MolGraph) -> Vector{Bool}

Return a vector of size ``n`` representing whether 1 to ``n``th bonds of
the given molecule belong to a ring or not.
"""
function is_edge_in_ring(mol::MolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :e_is_in_ring) && return get_cache(mol, :e_is_in_ring)
    return .!isempty.(edge_which_ring(mol))
end



"""
    bmscaffold(mol::MolGraph) -> Vector{Int}

Return Bemis-Murcko scaffold as a vector of scaffold vertices.
"""
function bmscaffold(g::SimpleGraph)
    (length(connected_components(g)) == 1
        || error("bmscaffold for a disconnected graph cannot be determined"))
    degree_ = degree(g)
    while true
        terms = findall(degree_ .== 1)
        isempty(terms) && break
        for t in terms
            degree_[t] = 0
            for nbr in neighbors(g, t)
                degree_[nbr] == 0 && continue
                degree_[nbr] -= 1
            end
        end
    end
    return findall(degree_ .!= 0)
end
bmscaffold(mol::SimpleMolGraph) = bmscaffold(mol.graph)


"""
    scaffold_fragments(g::SimpleGraph) -> Vector{Int}

Return scaffold fragments. Intended for a MMP derivative.
"""
function scaffold_fragments(g::SimpleGraph)
    bms = bmscaffold(g)
    isempty(bms) && return scaffolds  # no scaffold
    queue = [bms]  # BFS
    scaffolds = Vector{Int}[bms]
    dag = SimpleDiGraph{Int}(1)  # fragmentation relationship (parent -> child)
    while !isempty(queue)
        # generate new fragments
        scaffold = pop!(queue)
        scf, scfvmap = induced_subgraph(g, scaffold)
        parent =findfirst(x -> x == scaffold, scaffolds)
        for fr in fused_rings(scf)
            # enumerate connections to the fused ring
            nbrs = setdiff(union([neighbors(scf, v) for v in fr]...), fr)
            length(nbrs) == 1 || continue
            # remove a fused ring
            subg, vmap = induced_subgraph(g, setdiff(scaffold, scfvmap[fr]))
            fragment = vmap[bmscaffold(subg)]
            if fragment in scaffolds  # duplicate, skip the branch
                dst = findfirst(x -> x == fragment, scaffolds)
                add_edge!(dag, parent, dst)
                continue
            end
            pushfirst!(queue, fragment)
            push!(scaffolds, fragment)
            add_vertex!(dag)
            add_edge!(dag, parent, length(scaffolds))
            if length(scaffolds) == 400
                @info "max fragment count reached"
                break
            end
        end
    end
    return scaffolds, dag
end



# Primary properties

"""
    atom_symbol(mol::MolGraph) -> Vector{Symbol}

Return a vector of size ``n`` representing atom symbols of 1 to ``n``th atoms of
the given molecule.
"""
function atom_symbol(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_symbol) && return get_cache(mol, :v_symbol)
    return [get_prop(mol, i, :symbol) for i in vertices(mol)]
end

atom_symbol!(mol::SimpleMolGraph
    ) = set_cache!(mol, :v_symbol, [get_prop(mol, i, :symbol) for i in vertices(mol)])

"""
    charge(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing atom charges of 1 to ``n``th atoms of
the given molecule.
"""
function charge(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_charge) && return get_cache(mol, :v_charge)
    return [get_prop(mol, i, :charge) for i in vertices(mol)]
end

# kekulize! or charge standardization (e.g. polarize!) would be reset
charge!(mol::SimpleMolGraph
    ) = set_cache!(mol, :v_charge, [get_prop(mol, i, :charge) for i in vertices(mol)])


"""
    multiplicity(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing atom multiplicities of 1 to ``n``th atoms of
the given molecule (1: non-radical, 2: radical, 3: biradical).
"""
function multiplicity(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_multiplicity) && return get_cache(mol, :v_multiplicity)
    return [get_prop(mol, i, :multiplicity) for i in vertices(mol)]
end


"""
    bond_order(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing bond order of 1 to ``n``th bonds of
the given molecule.
"""
function bond_order(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :e_order) && return get_cache(mol, :e_order)
    return [get_prop(mol, e, :order) for e in edges(mol)]
end

# mass -> src/mass.jl
# coords -> src/coords.jl


# Secondary properties (descriptors)

# Valence

"""
    lone_pair(mol::MolGraph) -> Vector{Union{Int,Nothing}}

Return a vector of size ``n`` representing the number of lone pairs of
1 to ``n``th atoms of the given molecule.

The number of lone pair in inorganic atoms would be `nothing`.
The result can take negative value if the atom has empty shells (e.g. B).
"""
function lone_pair(symbol_arr, charge_arr)
    arr = Vector{Union{Int,Nothing}}(nothing, length(symbol_arr))
    for i in 1:length(symbol_arr)
        haskey(LONEPAIR_COUNT, symbol_arr[i]) || continue
        arr[i] = LONEPAIR_COUNT[symbol_arr[i]] - charge_arr[i]
    end
    return arr
end

function lone_pair(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_lone_pair) && return get_cache(mol, :v_lone_pair)
    return lone_pair(atom_symbol(mol), charge(mol))
end

lone_pair!(mol::SimpleMolGraph
    ) = set_cache!(mol, :v_lone_pair, lone_pair(atom_symbol(mol), charge(mol)))


"""
    apparent_valence(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of total bond order
incident to 1 to ``n``th atoms of the given molecule.
"""
function apparent_valence(g, order_arr)
    arr = fill(zero(Int), nv(g))
    er = Dict(e => i for (i, e) in enumerate(edges(g)))
    for e in edges(g)
        arr[src(e)] += order_arr[er[e]]
        arr[dst(e)] += order_arr[er[e]]
    end
    return arr
end

function apparent_valence(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_apparent_valence) && return get_cache(mol, :v_apparent_valence)
    return apparent_valence(mol.graph, bond_order(mol))
end

apparent_valence!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_apparent_valence, apparent_valence(mol.graph, bond_order(mol)))


"""
    valence(mol::SimpleMolGraph) -> Vector{Union{Int,Nothing}}

Return a vector of size ``n`` representing the intrinsic valence of
1 to ``n``th atoms of the given molecule.

The number of implicit hydrogens would be calculated based on the valence.
The valence value in inorganic atoms would be `nothing`.
This property corresponds to SMARTS `v` query.
"""
function valence(symbol_arr, lone_pair_arr)
    arr = Vector{Union{Int,Nothing}}(nothing, length(symbol_arr))
    for i in 1:length(symbol_arr)
        lone_pair_arr[i] === nothing && continue
        arr[i] = symbol_arr[i] === :H ? 1 : 4 - abs(lone_pair_arr[i])
    end
    return arr
end

function valence(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_valence) && return get_cache(mol, :v_valence)
    return valence(atom_symbol(mol), lone_pair(mol))
end

valence!(mol::SimpleMolGraph
    ) = set_cache!(mol, :v_valence, valence(atom_symbol(mol), lone_pair(mol)))


"""
    explicit_hydrogens(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of explicit hydrogens
connected to 1 to ``n``th atoms of the given molecule.

"Explicit" means hydrogens are explicitly represented as graph nodes.
"""
function explicit_hydrogens(g, symbol_arr)
    arr = fill(zero(Int), nv(g))
    for i in vertices(g)
        symbol_arr[i] === :H || continue
        for nbr in neighbors(g, i)
            arr[nbr] += 1
        end
    end
    return arr
end

function explicit_hydrogens(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_explicit_hydrogens) && return get_cache(mol, :v_explicit_hydrogens)
    return explicit_hydrogens(mol.graph, atom_symbol(mol))
end

explicit_hydrogens!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_explicit_hydrogens, explicit_hydrogens(mol.graph, atom_symbol(mol)))


"""
    implicit_hydrogens(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of implicit hydrogens
connected to 1 to ``n``th atoms of the given molecule.

"Implicit" means hydrogens are not represented as graph nodes,
but it can be infered from the intrinsic valence of typical organic atoms.
"""
function implicit_hydrogens(valence_arr, app_valence_arr)
    arr = fill(zero(Int), length(valence_arr))
    for i in 1:length(valence_arr)
        valence_arr[i] === nothing && continue
        arr[i] = max(0, valence_arr[i] - app_valence_arr[i])
    end
    return arr
end

function implicit_hydrogens(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_implicit_hydrogens) && return get_cache(mol, :v_implicit_hydrogens)
    return implicit_hydrogens(valence(mol), apparent_valence(mol))
end

implicit_hydrogens!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_implicit_hydrogens, implicit_hydrogens(valence(mol), apparent_valence(mol)))


"""
    heavyatoms(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of non-hydrogen atoms
connected to 1 to ``n``th atoms of the given molecule.
"""
function heavy_atoms(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_heavy_atoms) && return get_cache(mol, :v_heavy_atoms)
    return degree(mol) - explicit_hydrogens(mol)
end

heavy_atoms!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_heavy_atoms, degree(mol) - explicit_hydrogens(mol))



"""
    total_hydrogens(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of total hydrogens
(implicit and explicit) connected to 1 to ``n``th atoms of the given molecule.

This property corresponds to SMARTS `H` query.
"""
function total_hydrogens(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_total_hydrogens) && return get_cache(mol, :v_total_hydrogens)
    return explicit_hydrogens(mol) + implicit_hydrogens(mol)
end

total_hydrogens!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_total_hydrogens, explicit_hydrogens(mol) + implicit_hydrogens(mol))


"""
    connectivity(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of total atoms
(implicit and explicit) connected to 1 to ``n``th atoms of the given molecule.

This property corresponds to SMARTS `X` query.
"""
function connectivity(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_connectivity) && return get_cache(mol, :v_connectivity)
    return degree(mol) + implicit_hydrogens(mol)
end

connectivity!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_connectivity, degree(mol) + implicit_hydrogens(mol))



# Hydrogen bond donor/acceptor

function is_hydrogen_acceptor(symbol_arr, lone_pair_arr)
    ac = (sym, lp) -> lp === nothing ? false : sym in (:N, :O, :F) && lp > 0
    return ac.(symbol_arr, lone_pair_arr)
end

function is_hydrogen_acceptor(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_is_hydrogen_acceptor) && return get_cache(mol, :v_is_hydrogen_acceptor)
    return is_hydrogen_acceptor(atom_symbol(mol), lone_pair(mol))
end

is_hydrogen_acceptor!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_is_hydrogen_acceptor, is_hydrogen_acceptor(atom_symbol(mol), lone_pair(mol)))


"""
    hacceptorcount(mol::SimpleMolGraph) -> Int

Return the total number of hydrogen bond acceptors (N, O and F).
"""
hydrogen_acceptor_count(mol::SimpleMolGraph) = reduce(+, is_hydrogen_acceptor(mol); init=0)


function is_hydrogen_donor(symbol_arr, total_hydrogens_arr)
    dc = (sym, h) -> sym in (:N, :O) && h > 0
    return dc.(symbol_arr, total_hydrogens_arr)
end

function is_hydrogen_donor(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_is_hydrogen_donor) && return get_cache(mol, :v_is_hydrogen_donor)
    return is_hydrogen_donor(atom_symbol(mol), total_hydrogens(mol))
end

is_hydrogen_donor!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_is_hydrogen_donor, is_hydrogen_donor(atom_symbol(mol), total_hydrogens(mol)))


"""
    hdonorcount(mol::SimpleMolGraph) -> Int

Return the total number of hydrogen bond donors (O and N attached to hydrogens).
"""
hydrogen_donor_count(mol::SimpleMolGraph) = reduce(+, is_hydrogen_donor(mol); init=0)




# Rotatable bonds

"""
    isrotatable(mol::SimpleMolGraph)

Return a vector of size ``n`` representing whether 1 to ``n``th bonds
of the given molecule are rotatable or not.
"""
function is_rotatable(edge_list, degree_arr, edge_in_ring_arr, order_arr)
    arr = Vector{Bool}(undef, length(edge_list))
    for (i, e) in enumerate(edge_list)
        arr[i] = (!edge_in_ring_arr[i] && order_arr[i] == 1
            && degree_arr[src(e)] != 1 && degree_arr[dst(e)] != 1)
    end
    return arr
end

function is_rotatable(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :e_is_rotatable) && return get_cache(mol, :e_is_rotatable)
    return is_rotatable(edges(mol), degree(mol), is_edge_in_ring(mol), bond_order(mol))
end

is_rotatable!(mol::SimpleMolGraph) = set_cache!(
    mol, :e_is_rotatable,
    is_rotatable(edges(mol), degree(mol), is_edge_in_ring(mol), bond_order(mol)))


"""
    rotatablecount(mol::SimpleMolGraph) -> Int

Return the total number of rotatable bonds.
"""
rotatable_count(mol::SimpleMolGraph) = reduce(+, is_rotatable(mol); init=0)




# Composition

"""
    atomcounter(mol::SimpleMolGraph) -> Dict{Symbol,Int}

Count the number of atoms and return symbol => count dict.
"""
function atom_counter(mol::SimpleMolGraph)
    counter = Dict{Symbol,Int}()
    for sym in atom_symbol(mol)
        if !haskey(counter, sym)
            counter[sym] = 0
        end
        counter[sym] += 1
    end
    hcnt = reduce(+, implicit_hydrogens(mol); init=0)
    if hcnt > 0
        if !haskey(counter, :H)
            counter[:H] = 0
        end
        counter[:H] += hcnt
    end
    return counter
end


"""
    heavy_atom_count(mol::SimpleMolGraph) -> Int

Return the total number of non-hydrogen atoms.
"""
function heavy_atom_count(mol::SimpleMolGraph)
    counter = atom_counter(mol)
    delete!(counter, :H)
    return reduce(+, values(counter); init=0)
end


function write_formula(counter::Dict{Symbol,Int})
    strs = []
    if haskey(counter, :C)
        push!(strs, "C")
        c = pop!(counter, :C)
        c > 1 && push!(strs, string(c))
        if haskey(counter, :H)
            push!(strs, "H")
            h = pop!(counter, :H)
            h > 1 && push!(strs, string(h))
        end
    end
    for sym in sort(collect(keys(counter)))
        push!(strs, string(sym))
        cnt = counter[sym]
        cnt > 1 && push!(strs, string(cnt))
    end
    return join(strs)
end


"""
    molecular_formula(mol::MolGraph) -> String

Return the molecular formula in Hill system.
"""
function molecular_formula(mol::SimpleMolGraph)
    counter = atom_counter(mol)
    return write_formula(counter)
end


"""
    empirical_formula(mol::MolGraph) -> String

Return the empirical formula in Hill system.
"""
function empirical_formula(mol::SimpleMolGraph)
    counter = atom_counter(mol)
    if length(counter) > 1
        dv = reduce(gcd, values(counter))
        for k in keys(counter)
            counter[k] = div(counter[k], dv)
        end
    end
    return write_formula(counter)
end




# Hybridization

"""
    pi_electron(mol::SimpleMolGraph) -> Vector{Int}

Returns a vector of size ``n`` representing the number of ``\\pi`` electrons
of 1 to ``n``th atoms of the given molecule.

The counting of ``\\pi`` electrons is based on the following rules.

- sp3 (lone pair + connectivity == 4) or not organic heavy atom -> 0
- sp2 (lone pair + connectivity == 3) -> 1
- sp (lone pair + connectivity == 2) -> 2

Electron delocalization is not considered (e.g. mesomeric effect in carbonyl group).
"""
function pi_electron(lone_pair_arr, connectivity_arr)
    arr = fill(zero(Int), length(lone_pair_arr))
    @inbounds for i in 1:length(lone_pair_arr)
        if !isnothing(lone_pair_arr[i])
            orbitals = connectivity_arr[i] + lone_pair_arr[i]
            if orbitals == 3
                arr[i] = 1
            elseif orbitals == 2
                arr[i] = 2
            end
        end
    end
    return arr
end

function pi_electron(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_pi_electron) && return get_cache(mol, :v_pi_electron)
    return pi_electron(lone_pair(mol), connectivity(mol))
end

pi_electron!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_pi_electron, pi_electron(lone_pair(mol), connectivity(mol)))


"""
    pi_delocalized(mol::SimpleMolGraph) -> Vector{Int}

Returns a vector of size ``n`` representing the number of ``\\pi`` electrons
of 1 to ``n``th atoms of the given molecule.

The following delocalization rules are added to `pi_electron`.

- N (except for ammonium) or O adjacent to sp2 atom -> 2
"""
function pi_delocalized(g, symbol_arr, charge_arr, pi_arr)
    arr = copy(pi_arr)
    @inbounds for i in 1:length(pi_arr)
        if symbol_arr[i] in (:O, :N) && pi_arr[i] == 0 && charge_arr[i] <= 0
            for nbr in neighbors(g, i)
                if pi_arr[nbr] > 0
                    arr[i] = 2
                    break
                end
            end
        end
    end
    return arr
end

function pi_delocalized(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_pi_delocalized) && return get_cache(mol, :v_pi_delocalized)
    return pi_delocalized(mol.graph, atomsymbol(mol), charge(mol), pi_electron(mol))
end

pi_delocalized!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_pi_delocalized,
    pi_delocalized(mol.graph, atomsymbol(mol), charge(mol), pi_electron(mol)))

"""
    hybridization(mol::SimpleMolGraph) -> Vector{Int}

Returns a vector of size ``n`` representing the orbital hybridization symbols
(`:sp3`, `:sp2`, `:sp` or `:none`) of 1 to ``n``th atoms of the given molecule.

The hybridization value in inorganic atoms and non-typical organic atoms will be `:none`
(e.g. s, sp3d and sp3d2 orbitals). The categorization rule is the same as `pi_electron`,
but this method returns orbital hybridization symbols.
Electron delocalization is not considered (e.g. mesomeric effect in carbonyl group).
"""
function hybridization(lone_pair_arr, connectivity_arr)
    arr = fill(:none, length(lone_pair_arr))
    @inbounds for i in 1:length(lone_pair_arr)
        if !isnothing(lone_pair_arr[i])
            orbitals = connectivity_arr[i] + lone_pair_arr[i]
            if orbitals == 4
                arr[i] = :sp3
            elseif orbitals == 3
                arr[i] = :sp2
            elseif orbitals == 2
                arr[i] = :sp
            end
        end
    end
    return arr
end

function hybridization(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_hybridization) && return get_cache(mol, :v_hybridization)
    return hybridization(lone_pair(mol), connectivity(mol))
end

hybridization!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_hybridization, hybridization(lone_pair(mol), connectivity(mol)))

"""
    hybridization_delocalized(mol::SimpleMolGraph) -> Vector{Int}

Returns a vector of size ``n`` representing the number of ``\\pi`` electrons
of 1 to ``n``th atoms of the given molecule.

The following delocalization rules are added to `hybridization`.

- N (except for ammonium) or O adjacent to sp2 atom -> :sp2
"""
function hybridization_delocalized(g, symbol_arr, charge_arr, hyb_arr)
    arr = copy(hyb_arr)
    @inbounds for i in 1:length(hyb_arr)
        if symbol_arr[i] in (:O, :N) && hyb_arr[i] === :sp3 && charge_arr[i] <= 0
            for nbr in neighbors(g, i)
                if hyb_arr[nbr] in (:sp, :sp2)
                    arr[i] = :sp2
                    break
                end
            end
        end
    end
    return arr
end

function hybridization_delocalized(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_hybridization_delocalized) && return get_cache(mol, :v_hybridization_delocalized)
    return hybridization_delocalized(mol.graph, atomsymbol(mol), charge(mol), hybridization(mol))
end

hybridization_delocalized!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_hybridization_delocalized,
    hybridization_delocalized(mol.graph, atomsymbol(mol), charge(mol), hybridization(mol)))

# Aromaticity

function is_ring_aromatic(g, degree_arr, sssr_, symbol_arr, smiles_isaromatic_arr, lone_pair_arr, app_valence_arr)
    pie_ = app_valence_arr - degree_arr
    carbonyl_o = findall(
        (symbol_arr .=== :O) .* (degree_arr .== 1) .* (pie_ .== 1))
    carbonyl_c = falses(nv(g))
    for o in carbonyl_o
        c = neighbors(g, o)[1]
        if symbol_arr[c] === :C
            carbonyl_c[c] = true
        end
    end
    arr = falses(length(sssr_))
    for (i, ring) in enumerate(sssr_)
        if all(smiles_isaromatic_arr[ring])  # SMILES aromatic atom
            arr[i] = true
            continue
        end
        cnt = 0
        for r in ring
            carbonyl_c[r] && continue
            if pie_[r] == 1  # pi electron x1
                cnt += 1
                continue
            elseif lone_pair_arr[r] !== nothing
                if lone_pair_arr[r] > 0  # lone pair (pi electron x2)
                    cnt += 2  
                    continue
                elseif lone_pair_arr[r] < 0  # :B (pi electron x0)
                    continue  
                end
                # lonepair == 0
            end
            cnt = 0  # others are non-conjugated, cannot be aromatic
            break
        end
        if cnt % 4 == 2  # Huckel rule check
            arr[i] = true
        end
    end
    return arr
end

function is_ring_aromatic(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :is_ring_aromatic) && return get_cache(mol, :is_ring_aromatic)
    arom_arr = (hasfield(vproptype(mol), :isaromatic)
        ? [get_prop(mol, i, :isaromatic) for i in vertices(mol)] : falses(nv(mol)))
    return is_ring_aromatic(
        mol.graph, degree(mol), sssr(mol), atom_symbol(mol), arom_arr, lone_pair(mol), apparent_valence(mol))
end

function is_ring_aromatic!(mol::SimpleMolGraph)
    arom_arr = (hasfield(vproptype(mol), :isaromatic)
        ? [get_prop(mol, i, :isaromatic) for i in vertices(mol)] : falses(nv(mol)))
    set_cache!(mol, :is_ring_aromatic, is_ring_aromatic(
        mol.graph, degree(mol), sssr(mol), atom_symbol(mol), arom_arr, lone_pair(mol), apparent_valence(mol)))
end


"""
    is_aromatic(mol::SimpleMolGraph) -> Vector{Bool}

Returns a vector of size ``n`` representing whether 1 to ``n``th atoms
of the given molecule belong to an aromatic ring or not.

Some kind of aromaticity resulting from long conjugated chains and charge
delocalization may be unrecognizable. Also, non-classical aromaticity
such as Moebius aromaticity is not considered.
"""
function is_aromatic(g, sssr_, is_ring_arom)
    arr = falses(nv(g))
    for ring in sssr_[findall(is_ring_arom)]
        arr[ring] .= true
    end
    return arr
end


function is_aromatic(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_is_aromatic) && return get_cache(mol, :v_is_aromatic)
    return is_aromatic(mol.graph, sssr(mol), is_ring_aromatic(mol))
end

is_aromatic!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_is_aromatic, is_aromatic(mol.graph, sssr(mol), is_ring_aromatic(mol)))


"""
    is_edge_aromatic(mol::SimpleMolGraph) -> Vector{Bool}

Returns a vector of size ``n`` representing whether 1 to ``n``th bonds
of the given molecule belong to an aromatic ring or not.

See [`is_aromatic`](@ref).
"""
function is_edge_aromatic(g, sssr_, is_ring_arom)
    arr = falses(ne(g))
    er = Dict(e => i for (i, e) in enumerate(edges(g)))
    for ring in sssr_[findall(is_ring_arom)]
        for i in 1:(length(ring) - 1)
            arr[er[u_edge(g, ring[i], ring[i + 1])]] = true
        end
        arr[er[u_edge(g, ring[1], ring[end])]] = true
    end
    return arr
end

function is_edge_aromatic(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :e_is_aromatic) && return get_cache(mol, :e_is_aromatic)
    return is_edge_aromatic(mol.graph, sssr(mol), is_ring_aromatic(mol))
end

is_edge_aromatic!(mol::SimpleMolGraph) = set_cache!(
    mol, :e_is_aromatic, is_edge_aromatic(mol.graph, sssr(mol), is_ring_aromatic(mol)))



# deprecated function names

nodedegree = degree
sssrmembership = which_ring
sssrbondmembership = edge_which_ring
fusedrings = fused_rings
fusedringmembership = which_fused_ring
smallestsssr = smallest_ring
sssrcount = ring_count
isringatom = is_in_ring
isringbond = is_edge_in_ring
atomsymbol(mol::MolGraph) = atom_symbol(mol)
bondorder = bond_order
lonepair = lone_pair
heavyatomconnected = heavy_atoms
explicithconnected = explicit_hydrogens
implicithconnected = implicit_hydrogens
hydrogenconnected = total_hydrogens
ishdonor = is_hydrogen_donor
hdonorcount = hydrogen_donor_count
ishacceptor = is_hydrogen_acceptor
hacceptorcount = hydrogen_acceptor_count
isrotatable = is_rotatable
rotatablecount = rotatable_count
atomcounter = atom_counter
heavyatomcount = heavy_atom_count
molecularformula = molecular_formula
empiricalformula = empirical_formula
pielectron = pi_electron
isaromaticring = is_ring_aromatic
isaromatic = is_aromatic
isaromaticbond = is_edge_aromatic