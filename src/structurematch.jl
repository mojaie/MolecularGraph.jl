#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# Return `string(descriptors[key][i])`, or `nothing` if `key` is not in `descriptors`.
# Unrolled at compile time so that each descriptor array keeps its concrete type.
@generated function descriptor_string(
        descriptors::NamedTuple{K}, key::Symbol, i::Integer) where K
    ex = :(nothing)
    for k in reverse(K)
        ex = :(key === $(QuoteNode(k)) ? string(descriptors.$k[i]) : $ex)
    end
    return ex
end


"""
    vmatchgen(mol1::MolGraph, mol2::MolGraph) -> Function
    vmatchgen(mol1::MolGraph{T1,V1,E1}, mol2::MolGraph{T2,V2,E2}
        ) where {T1,T2,V1,V2<:QueryTree,E1,E2} -> Function
    vmatchgen(mol1::MolGraph{T1,V1,E1}, mol2::MolGraph{T2,V2,E2}
        ) where {T1,T2,V1<:QueryTree,V2<:QueryTree,E1,E2}

Return a default vertex attribute matching function for graph isomorphism algorithms.
"""
function vmatchgen(mol1::SimpleMolGraph, mol2::SimpleMolGraph)
    sym1 = atom_symbol(mol1)
    sym2 = atom_symbol(mol2)
    geo1 = hybridization(mol1)
    geo2 = hybridization(mol2)
    return (v1, v2) -> sym1[v1] == sym2[v2] && geo1[v1] === geo2[v2]
end

function vmatchvecgen(mol::SimpleMolGraph)
    sym = atom_symbol(mol)  # atom number 1-118
    geo = hybridization(mol)
    conv = Dict(:none => 0, :sp => 1, :sp2 => 2, :sp3 => 3)
    return function (i)
        return conv[geo[i]] * 120 + atom_number(sym[i])
    end
end

function vmatchgen(mol1::SimpleMolGraph, mol2::T) where T <: QueryMolGraph
    descriptors = (  # precalculated descriptors
        symbol = atom_symbol(mol1),
        isaromatic = is_aromatic(mol1),
        charge = atom_charge(mol1),
        isotope = isotope(mol1),
        connectivity = connectivity(mol1),
        degree = degree(mol1),
        valence = valence(mol1),
        total_hydrogens = total_hydrogens(mol1),
        smallest_ring = smallest_ring(mol1),
        ring_count = ring_count(mol1),
    )
    recursive = Dict{String,T}()  # cache recursive mols
    matches = Dict{Int,Dict{Int,Bool}}()  # cache matches
    return function (v1, v2)
        haskey(matches, v1) && haskey(matches[v1], v2) && return matches[v1][v2]
        qtree = mol2[v2]
        nodes = QueryNode[
            qand(),
            qeq(:symbol, descriptor_string(descriptors, :symbol, v1)),
            qeq(:isaromatic, descriptor_string(descriptors, :isaromatic, v1)),
            qeq(:charge, descriptor_string(descriptors, :charge, v1)),
            qeq(:isotope, descriptor_string(descriptors, :isotope, v1)),
            qeq(:connectivity, descriptor_string(descriptors, :connectivity, v1)),
            qeq(:degree, descriptor_string(descriptors, :degree, v1)),
            qeq(:valence, descriptor_string(descriptors, :valence, v1)),
            qeq(:total_hydrogens, descriptor_string(descriptors, :total_hydrogens, v1)),
            qeq(:smallest_ring, descriptor_string(descriptors, :smallest_ring, v1)),
            qeq(:ring_count, descriptor_string(descriptors, :ring_count, v1))
        ]
        edges = [(1, i) for i in 2:length(nodes)]
        # Recursive queries depend on the environment of the atom in mol1, which cannot be
        # inferred from the descriptors above. Evaluate them by substructure match and
        # add the results as literals (e.g. $(NC=O) or !$(NC=O)).
        for p in unique(get(querypropmap(qtree), :recursive, QueryNode[]))
            haskey(recursive, p.value) || (recursive[p.value] = smartstomol(p.value))
            if has_substruct_match(mol1, recursive[p.value], mandatory=Dict(v1 => 1))
                push!(nodes, p)
                push!(edges, (1, length(nodes)))
            else
                push!(nodes, qnot())
                push!(edges, (1, length(nodes)))
                push!(nodes, p)
                push!(edges, (length(nodes) - 1, length(nodes)))
            end
        end
        res = issubset(QueryAtom(edges, nodes), qtree)
        haskey(matches, v1) || (matches[v1] = Dict{Int,Bool}())
        matches[v1][v2] = res
        return res
    end
end

function vmatchgen(mol1::QueryMolGraph, mol2::QueryMolGraph)
    recursive = Dict{String,Dict{String,Bool}}()  # cache recursive queries
    matches = Dict{Int,Dict{Int,Bool}}()  # cache matches
    return function (v1, v2)
        haskey(matches, v1) && haskey(matches[v1], v2) && return matches[v1][v2]
        res = issubset(mol1[v1], mol2[v2]; recursive_cache=recursive)
        haskey(matches, v1) || (matches[v1] = Dict{Int,Bool}())
        matches[v1][v2] = res
        return res
    end
end


"""
    ematchgen(mol1::MolGraph, mol2::MolGraph) -> Function
    ematchgen(mol::MolGraph{T1,V1,E1}, qmol::MolGraph{T2,V2,E2}
        ) where {T1,T2,V1,V2,E1,E2<:QueryTree}
    ematchgen(mol1::MolGraph{T1,V1,E1}, mol2::MolGraph{T2,V2,E2}
        ) where {T1,T2,V1,V2,E1<:QueryTree,E2<:QueryTree}

Return a default edge attribute matching function for graph isomorphism algorithms.
"""
function ematchgen(mol1::SimpleMolGraph, mol2::SimpleMolGraph)
    return (e1, e2) -> true
end

function ematchvecgen(mol::SimpleMolGraph)
    return (u, v) -> 0
end

function ematchgen(mol1::SimpleMolGraph, mol2::QueryMolGraph)
    descriptors = (  # precalculated descriptors
        order = bond_order(mol1),
        is_in_ring = is_edge_in_ring(mol1),
        isaromatic = is_edge_aromatic(mol1)
    )
    matches = Dict{Edge{Int},Dict{Edge{Int},Bool}}()  # cache matches
    ernk1 = edge_rank(mol1)
    return function (e1, e2)
        haskey(matches, e1) && haskey(matches[e1], e2) && return matches[e1][e2]
        e1q = QueryBond(
            [(1, 2), (1, 3), (1, 4)],
            [
                qand(),
                qeq(:order, descriptor_string(descriptors, :order, ernk1[e1])),
                qeq(:is_in_ring, descriptor_string(descriptors, :is_in_ring, ernk1[e1])),
                qeq(:isaromatic, descriptor_string(descriptors, :isaromatic, ernk1[e1]))
            ]
        )
        res = issubset(e1q, mol2[e2])
        haskey(matches, e1) || (matches[e1] = Dict{Edge{Int},Bool}())
        matches[e1][e2] = res
        return res
    end
end

function ematchgen(mol1::QueryMolGraph, mol2::QueryMolGraph)
    matches = Dict{Edge{Int},Dict{Edge{Int},Bool}}()  # cache matches
    return function (e1, e2)
        haskey(matches, e1) && haskey(matches[e1], e2) && return matches[e1][e2]
        res = issubset(mol1[e1], mol2[e2])
        haskey(matches, e1) || (matches[e1] = Dict{Edge{Int},Bool}())
        matches[e1][e2] = res
        return res
    end
end


function exact_topology_prefilter(mol1::SimpleMolGraph, mol2::SimpleMolGraph)
    nv(mol1) == nv(mol2) || return false
    ne(mol1) == ne(mol2) || return false
    if has_descriptor(mol1, Val(:sssr)) && has_descriptor(mol2, Val(:sssr))
        r1 = get_descriptor(mol1, Val(:sssr))
        r2 = get_descriptor(mol2, Val(:sssr))
        length(r1) == length(r2) || return false
    end
    return true
end

function topology_prefilter(mol1::SimpleMolGraph, mol2::SimpleMolGraph)
    nv(mol1) >= nv(mol2) || return false
    ne(mol1) >= ne(mol2) || return false
    if has_descriptor(mol1, Val(:sssr)) && has_descriptor(mol2, Val(:sssr))
        r1 = get_descriptor(mol1, Val(:sssr))
        r2 = get_descriptor(mol2, Val(:sssr))
        length(r1) >= length(r2) ||  return false
    end
    return true
end


function node_stereo_match(
        revmap::Dict{T,T}, # mol2 => mol1
        mol1::SimpleMolGraph{T}, mol2::SimpleMolGraph{T}) where T
    for (i, conf) in mol2[:stereocenter]
        f, s, t, dr = conf
        revmap[i] in keys(mol1[:stereocenter]) || return false
        isc = isclockwise(mol1[:stereocenter][revmap[i]], revmap[f], revmap[s], revmap[t])
        isc === dr || return false
    end
    for (e, conf) in mol2[:stereobond]
        f2, s2, iscis2 = conf
        e1 = u_edge(T, revmap[e.src], revmap[e.dst])
        e1 in keys(mol1[:stereobond]) || return false
        f1, s1, iscis1 = mol1[:stereobond][e1]
        fmatch = f1 in [revmap[f2], revmap[s2]]
        smatch = s1 in [revmap[f2], revmap[s2]]
        (fmatch === smatch) === (iscis1 === iscis2) || return false
    end
    return true
end

function node_stereo_match(
        revmap::Dict{T,T}, # mol2 => mol1
        mol1::SimpleMolGraph{T}, mol2::QueryMolGraph{T,V,E}) where {T,V,E}
    for (i, conf) in mol2[:stereocenter]
        f, s, t, dr = conf
        if !haskey(mol1[:stereocenter], revmap[i])
            for prop in values(mol2[i].vprops)
                prop.key === :stereo || continue
                if prop.value in ("clockwiseuns", "anticlockwiseuns")
                    return true  # allow unspecified
                end
                break
            end
            return false
        end
        isc = isclockwise(mol1[:stereocenter][revmap[i]], revmap[f], revmap[s], revmap[t])
        isc === dr || return false
    end
    for (e, conf) in mol2[:stereobond]
        s2, d2, iscis2 = conf
        e1 = u_edge(T, revmap[e.src], revmap[e.dst])
        if !haskey(mol1[:stereobond], e1)
            srcp = values(mol2[u_edge(T, e.src, s2)].vprops)
            dstp = values(mol2[u_edge(T, e.dst, d2)].vprops)
            for prop in vcat(collect(srcp), collect(dstp))
                prop.key === :direction || continue
                if prop.value in ("upuns", "downuns")
                    return true  # allow unspecified
                end
                break
            end
            return false
        end
        f1, s1, iscis1 = mol1[:stereobond][e1]
        fmatch = f1 in [revmap[s2], revmap[d2]]
        smatch = s1 in [revmap[s2], revmap[d2]]
        (fmatch === smatch) === (iscis1 === iscis2) || return false
    end
    return true
end

function node_stereo_match(
        revmap::Dict{T,T}, # mol2 => mol1
        mol1::QueryMolGraph{T,V,E}, mol2::QueryMolGraph{T,V,E}) where {T,V,E}
    for (i, conf) in mol2[:stereocenter]
        f, s, t, dr = conf
        haskey(mol1[:stereocenter], revmap[i]) || return false
        # "or unspecified" is more general than clockwise/anticlockwise
        for prop in values(mol2[i].vprops)
            prop.key === :stereo || continue
            prop.value in ("clockwise", "anticlockwise") || continue
            for p1 in values(mol1[revmap[i]].vprops)
                p1.key === :stereo || continue
                p1.value in ("clockwiseuns", "anticlockwiseuns") && return false
                break
            end
            break
        end
        isc = isclockwise(mol1[:stereocenter][revmap[i]], revmap[f], revmap[s], revmap[t])
        isc === dr || return false
    end
    for (e, conf) in mol2[:stereobond]
        s2, d2, iscis2 = conf
        e1 = u_edge(T, revmap[e.src], revmap[e.dst])
        haskey(mol1[:stereobond], e1) || return false
        # "or unspecified" is more general than up/down
        srcp = values(mol2[u_edge(T, e.src, s2)].vprops)
        dstp = values(mol2[u_edge(T, e.dst, d2)].vprops)
        for prop in vcat(collect(srcp), collect(dstp))
            prop.key === :direction || continue
            prop.value in ("up", "down") || continue
            srcp1 = values(mol1[u_edge(T, revmap[e.src], revmap[s2])].vprops)
            dstp1 = values(mol1[u_edge(T, revmap[e.dst], revmap[d2])].vprops)
            for p1 in vcat(collect(srcp1), collect(dstp1))
                p1.key === :direction || continue
                p1.value in ("upuns", "downuns") && return false
                break
            end
            break
        end
        f1, s1, iscis1 = mol1[:stereobond][e1]
        fmatch = f1 in [revmap[s2], revmap[d2]]
        smatch = s1 in [revmap[s2], revmap[d2]]
        (fmatch === smatch) === (iscis1 === iscis2) || return false
    end
    return true
end

"""
    exact_matches(mol1, mol2; kwargs...) -> Iterator

Return a lazy iterator that generate node mappings between `mol` and `query` if these are exactly same.
See [`substruct_matches`](@ref) for available options.
"""
function exact_matches(mol1::SimpleMolGraph, mol2::SimpleMolGraph;
        vmatchgen=vmatchgen, ematchgen=ematchgen, stereo=false, kwargs...)
    # Note: InChI is better if you don't need mapping
    (nv(mol1) == 0 || nv(mol2) == 0) && return ()
    exact_topology_prefilter(mol1, mol2) || return ()
    mappings = isomorphisms(
        mol1.graph, mol2.graph,
        vmatch=vmatchgen(mol1, mol2), ematch=ematchgen(mol1, mol2); kwargs...)
    if stereo
        return Iterators.filter(mappings) do mapping
            return node_stereo_match(Dict(v => i for (i, v) in mapping), mol1, mol2)
        end
    else
        return mappings
    end
end


"""
    has_exact_match(mol1, mol2; kwargs...) -> Bool

Return whether `mol` and `query` have exactly the same structure.
See [`substruct_matches`](@ref) for available options.
"""
has_exact_match(mol1, mol2; kwargs...) = !isempty(exact_matches(mol1, mol2; kwargs...))


"""
    substruct_matches(mol1, mol2; kwargs...) -> Iterator

Return a lazy iterator that generate node mappings between `mol` and `query` if `mol` has `query` as a substructure.

# options

- `vmatch::Function`: a function for semantic atom attribute matching (default: `MolecularGraph.vmatch`)
- `ematch::Function`: a function for semantic bond attribute matching (default: `MolecularGraph.ematch`)
- `mandatory::Dict{Int,Int}`: mandatory node mapping (or edge mapping if matchtype=:edgeinduced)
- `forbidden::Dict{Int,Int}`: forbidden node mapping (or edge mapping if matchtype=:edgeinduced)
- `timeout::Union{Int,Nothing}`: if specified, abort vf2 calculation when the time reached and return empty iterator (default: 10 seconds).

"""
function substruct_matches(mol1::SimpleMolGraph, mol2::SimpleMolGraph;
        vmatchgen=vmatchgen, ematchgen=ematchgen, stereo=false, kwargs...)
    (nv(mol1) == 0 || nv(mol2) == 0) && return ()
    topology_prefilter(mol1, mol2) || return ()
    mappings = subgraph_monomorphisms(
        mol1.graph, mol2.graph,
        vmatch=vmatchgen(mol1, mol2), ematch=ematchgen(mol1, mol2); kwargs...)
    if stereo
        return Iterators.filter(mappings) do mapping
            return node_stereo_match(Dict(v => i for (i, v) in mapping), mol1, mol2)
        end
    else
        return mappings
    end
end


"""
    has_substruct_match(mol1, mol2; kwargs...) -> Bool

Return whether `mol` has `query` as a substructure.
See [`substruct_matches`](@ref) for available options.
"""
has_substruct_match(mol1, mol2; kwargs...
    ) = !isempty(substruct_matches(mol1, mol2; kwargs...))


"""
    node_substruct_matches(mol1, mol2; kwargs...) -> Iterator

Return a lazy iterator that generate node mappings between `mol` and `query` if `mol` has `query` as a substructure.
See [`substruct_matches`](@ref) for available options.
"""
function node_substruct_matches(mol1::SimpleMolGraph, mol2::SimpleMolGraph;
        vmatchgen=vmatchgen, ematchgen=ematchgen, stereo=false, kwargs...)
    (nv(mol1) == 0 || nv(mol2) == 0) && return ()
    topology_prefilter(mol1, mol2) || return ()
    mappings = nodesubgraph_isomorphisms(
        mol1.graph, mol2.graph,
        vmatch=vmatchgen(mol1, mol2), ematch=ematchgen(mol1, mol2); kwargs...)
    if stereo
        return Iterators.filter(mappings) do mapping
            return node_stereo_match(Dict(v => i for (i, v) in mapping), mol1, mol2)
        end
    else
        return mappings
    end
end


"""
    has_node_substruct_match(mol1, mol2; kwargs...) -> Bool

Return whether `mol` has `query` as a substructure.
See [`substruct_matches`](@ref) for available options.
"""
has_node_substruct_match(mol1, mol2; kwargs...
    ) = !isempty(node_substruct_matches(mol1, mol2; kwargs...))


"""
    edge_substruct_matches(mol1, mol2; kwargs...) -> Iterator

Return a lazy iterator that generate node mappings between `mol` and `query` if `mol` has `query` as a substructure.
See [`substruct_matches`](@ref) for available options.
"""
function edge_substruct_matches(mol1::SimpleMolGraph, mol2::SimpleMolGraph;
        vmatchgen=vmatchgen, ematchgen=ematchgen, kwargs...)
    (ne(mol1) == 0 || ne(mol2) == 0) && return ()
    topology_prefilter(mol1, mol2) || return ()
    return edgesubgraph_isomorphisms(
        mol1.graph, mol2.graph,
        vmatch=vmatchgen(mol1, mol2), ematch=ematchgen(mol1, mol2),
        ggmatch=vmatchgen(mol1, mol1), hhmatch=vmatchgen(mol2, mol2); kwargs...)
end


"""
    has_edge_substruct_match(mol1, mol2; kwargs...) -> Bool

Return whether `mol` has `query` as a substructure.
See [`substruct_matches`](@ref) for available options.
"""
has_edge_substruct_match(mol1, mol2; kwargs...
    ) = !isempty(edge_substruct_matches(mol1, mol2; kwargs...))



# MCS


function mcis_constraints(
        mol::SimpleMolGraph, vmatchvecgen=vmatchvecgen, ematchvecgen=ematchvecgen)
    return mcis_constraints(
        Val{:connection}(), mol.graph,
        vmatchvec=vmatchvecgen(mol), ematchvec=ematchvecgen(mol))
end

function mces_constraints(
        mol::SimpleMolGraph, vmatchvecgen=vmatchvecgen, ematchvecgen=ematchvecgen)
    return mces_constraints(
        Val{:connection}(), mol.graph,
        vmatchvec=vmatchvecgen(mol), ematchvec=ematchvecgen(mol))
end

function tdmcis_constraints(
        mol::SimpleMolGraph; vmatchvecgen=vmatchvecgen,
        ematchvecgen=ematchvecgen, kwargs...)
    return mcis_constraints(
        Val{:shortest}(), mol.graph,
        vmatchvec=vmatchvecgen(mol), ematchvec=ematchvecgen(mol); kwargs...)
end
function tdmces_constraints(
        mol::SimpleMolGraph; vmatchvecgen=vmatchvecgen,
        ematchvecgen=ematchvecgen, kwargs...)
    return mces_constraints(
        Val{:shortest}(), mol.graph,
        vmatchvec=vmatchvecgen(mol), ematchvec=ematchvecgen(mol); kwargs...)
end



"""
    disconnected_mcis(mol1::SimpleMolGraph, mol2::SimpleMolGraph; kwargs...) -> (Dict, Symbol)
    disconnected_mces(mol1::SimpleMolGraph, mol2::SimpleMolGraph; kwargs...) -> (Dict, Symbol)

Compute disconnected maximum common substructure (MCS) of mol1 and mol2.

## Keyword arguments

- timeout(Real): abort calculation and return suboptimal results if the execution
time has reached the given value (default=60, in seconds).
- targetsize(Int): abort calculation and return suboptimal result so far if the
given mcs size achieved.
"""
function disconnected_mcis(mol1::SimpleMolGraph, mol2::SimpleMolGraph; kwargs...)
    mol1 = mcis_constraints(mol1)
    mol2 = mcis_constraints(mol2)
    return maximum_common_subgraph(mol1, mol2; kwargs...)
end

function disconnected_mces(mol1::SimpleMolGraph, mol2::SimpleMolGraph; kwargs...)
    mol1 = mces_constraints(mol1)
    mol2 = mces_constraints(mol2)
    return maximum_common_subgraph(mol1, mol2; kwargs...)
end


"""
    connected_mcis(mol1::SimpleMolGraph, mol2::SimpleMolGraph; kwargs...) -> (Dict, Symbol)
    connected_mces(mol1::SimpleMolGraph, mol2::SimpleMolGraph; kwargs...) -> (Dict, Symbol)

Compute connected maximum common substructure (MCS) of mol1 and mol2.

## Keyword arguments

- timeout(Real): abort calculation and return suboptimal results if the execution
time has reached the given value (default=60, in seconds).
- targetsize(Int): abort calculation and return suboptimal result so far if the
given mcs size achieved.
"""
function connected_mcis(mol1::SimpleMolGraph, mol2::SimpleMolGraph; kwargs...)
    mol1 = mcis_constraints(mol1)
    mol2 = mcis_constraints(mol2)
    return maximum_common_subgraph(mol1, mol2, connected=true; kwargs...)
end

function connected_mces(mol1::SimpleMolGraph, mol2::SimpleMolGraph; kwargs...)
    mol1 = mces_constraints(mol1)
    mol2 = mces_constraints(mol2)
    return maximum_common_subgraph(mol1, mol2, connected=true; kwargs...)
end


"""
    tdmcis(mol1::SimpleMolGraph, mol2::SimpleMolGraph; kwargs...) -> (Dict, Symbol)
    tdmces(mol1::SimpleMolGraph, mol2::SimpleMolGraph; kwargs...) -> (Dict, Symbol)

Compute disconnected MCS of mol1 and mol2 with topological constraint (td-MCS).

## Keyword arguments

- diameter(Int): distance cutoff for topological constraint.
- tolerance(Int): distance mismatch tolerance for topological constraint.
- timeout(Real): abort calculation and return suboptimal results if the execution
time has reached the given value (default=60, in seconds).
- targetsize(Int): abort calculation and return suboptimal result so far if the
given mcs size achieved.

# References

1. Kawabata, T. (2011). Build-Up Algorithm for Atomic Correspondence between
Chemical Structures. Journal of Chemical Information and Modeling, 51(8),
1775–1787. https://doi.org/10.1021/ci2001023
1. https://www.jstage.jst.go.jp/article/ciqs/2017/0/2017_P4/_article/-char/en
"""
function tdmcis(mol1::SimpleMolGraph, mol2::SimpleMolGraph; kwargs...)
    mol1 = tdmcis_constraints(mol1; kwargs...)
    mol2 = tdmcis_constraints(mol2; kwargs...)
    return maximum_common_subgraph(mol1, mol2; kwargs...)
end

function tdmces(mol1::SimpleMolGraph, mol2::SimpleMolGraph; kwargs...)
    mol1 = tdmces_constraints(mol1; kwargs...)
    mol2 = tdmces_constraints(mol2; kwargs...)
    return maximum_common_subgraph(mol1, mol2; kwargs...)
end


"""
    nmap = emaptonmap(emap, mol, query)

Convert an edge-based mapping, of the form returned by [`edge_substruct_matches`](@ref), into
a map between nodes. Commonly, `nmap[i]` is a length-1 vector `[j]`, where `i=>j` is the mapping
from `nodeattr(query, i)` to `nodeattr(mol, j)`. In cases where the mapping is ambiguous,
`nmap[i]` may be multivalued.
"""
function emaptonmap(emap, mol::SimpleMolGraph, query::SimpleMolGraph)
    nmol, nq = nv(mol), nv(query)
    nq <= nmol || throw(ArgumentError("query must be a substructure of mol"))
    # Each node in the query edges can map to either of two nodes in mol
    qnodeoptions = [Tuple{Int,Int}[] for _ = 1:nq]
    for (edgemol, edgeq) in emap
        for nq in Tuple(edgeq)
            push!(qnodeoptions[nq], Tuple(edgemol))
        end
    end
    # For nodes connected to two or more other nodes, intersection results in a unique answer
    assignment = [intersect(nodeops...) for nodeops in qnodeoptions]
    # For the singly-connected nodes, assign them by eliminating ones already taken
    taken = falses(nmol)
    for a in assignment
        if length(a) == 1
            taken[a[1]] = true
        end
    end
    for a in assignment
        if length(a) > 1
            deleteat!(a, findall(taken[a]))
        end
    end
    return assignment
end
