#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

"""
    vmatchgen(mol1::MolGraph, mol2::MolGraph) -> Function
    vmatchgen(mol1::MolGraph{T1,V1,E1}, mol2::MolGraph{T2,V2,E2}
        ) where {T1,T2,V1,V2<:QueryTree,E1,E2} -> Function
    vmatchgen(mol1::MolGraph{T1,V1,E1}, mol2::MolGraph{T2,V2,E2}
        ) where {T1,T2,V1<:QueryTree,V2<:QueryTree,E1,E2}

Return a default vertex attribute matching function for graph isomorphism algorithms.
"""
function vmatchgen(mol1::MolGraph, mol2::MolGraph)
    sym1 = atom_symbol(mol1)
    sym2 = atom_symbol(mol2)
    geo1 = hybridization(mol1)
    geo2 = hybridization(mol2)
    return (v1, v2) -> sym1[v1] == sym2[v2] && geo1[v1] === geo2[v2]
end

function vmatchvecgen(mol)
    sym = atom_symbol(mol)  # atom number 1-118
    geo = hybridization(mol)
    conv = Dict(:none => 0, :sp => 1, :sp2 => 2, :sp3 => 3)
    return function (i)
        return conv[geo[i]] * 120 + atom_number(sym[i])
    end
end

function vmatchgen(mol1::MolGraph, mol2::T) where T <: QueryMolGraph
    descriptors = Dict(  # precalculated descriptors
        :symbol => atom_symbol(mol1),
        :isaromatic => is_aromatic(mol1),
        :charge => atom_charge(mol1),
        :mass => [atom_mass(props(mol1, i)) for i in vertices(mol1)],
        :connectivity => connectivity(mol1),
        :degree => degree(mol1),
        :valence => valence(mol1),
        :total_hydrogens => total_hydrogens(mol1),
        :smallest_ring => smallest_ring(mol1),
        :ring_count => ring_count(mol1),
    )
    recursive = Dict{String,T}()  # cache recursive mols
    matches = Dict{Int,Dict{Int,Bool}}()  # cache matches
    return function (v1, v2)
        haskey(matches, v1) && haskey(matches[v1], v2) && return matches[v1][v2]
        qtree = props(mol2, v2)
        qprop = union(QueryNode[], values(querypropmap(qtree))...)
        arr = falses(length(qprop))
        for (i, p) in enumerate(qprop)
            if p.key == :recursive
                haskey(recursive, p.value) || (recursive[p.value] = smartstomol(p.value))
                arr[i] = has_substruct_match(mol1, recursive[p.value], mandatory=Dict(v1 => 1))
            else
                arr[i] = string(descriptors[p.key][v1]) == p.value
            end
        end
        res = generate_queryfunc(qtree, qprop)(arr)
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
        res = issubset(props(mol1, v1), props(mol2, v2); recursive_cache=recursive)
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
function ematchgen(mol1::MolGraph, mol2::MolGraph)
    return (e1, e2) -> true
end

function ematchvecgen(mol)
    return (u, v) -> 0
end

function ematchgen(mol1::MolGraph, mol2::QueryMolGraph)
    descriptors = Dict(  # precalculated descriptors
        :order => bond_order(mol1),
        :is_in_ring => is_edge_in_ring(mol1),
        :isaromatic => is_edge_aromatic(mol1)
    )
    matches = Dict{Edge{Int},Dict{Edge{Int},Bool}}()  # cache matches
    return function (e1, e2)
        haskey(matches, e1) && haskey(matches[e1], e2) && return matches[e1][e2]
        qtree = props(mol2, e2)
        qprop = union(QueryNode[], values(querypropmap(qtree))...)
        arr = [string(descriptors[p.key][edge_rank(mol1, e1)]) == p.value for p in qprop]
        res = generate_queryfunc(qtree, qprop)(arr)
        haskey(matches, e1) || (matches[e1] = Dict{Edge{Int},Bool}())
        matches[e1][e2] = res
        return res
    end
end

function ematchgen(mol1::QueryMolGraph, mol2::QueryMolGraph)
    matches = Dict{Edge{Int},Dict{Edge{Int},Bool}}()  # cache matches
    return function (e1, e2)
        haskey(matches, e1) && haskey(matches[e1], e2) && return matches[e1][e2]
        res = issubset(props(mol1, e1), props(mol2, e2))
        haskey(matches, e1) || (matches[e1] = Dict{Edge{Int},Bool}())
        matches[e1][e2] = res
        return res
    end
end



function circuitrank(g::SimpleGraph)
    nv(g) == 0 && return 0
    return ne(g) - nv(g) + length(connected_components(g))
end


function exact_match_prefilter(mol1::SimpleMolGraph, mol2::SimpleMolGraph)
    nv(mol1) == nv(mol2) || return false
    ne(mol1) == ne(mol2) || return false
    circuitrank(mol1.graph) == circuitrank(mol2.graph) ||  return false
    return true
end


function substruct_match_prefilter(mol1::SimpleMolGraph, mol2::SimpleMolGraph)
    nv(mol1) >= nv(mol2) || return false
    ne(mol1) >= ne(mol2) || return false
    circuitrank(mol1.graph) >= circuitrank(mol2.graph) ||  return false
    return true
end


"""
    exact_matches(mol1, mol2; kwargs...) -> Iterator

Return a lazy iterator that generate node mappings between `mol` and `query` if these are exactly same.
See [`substruct_matches`](@ref) for available options.
"""
function exact_matches(mol1::SimpleMolGraph, mol2::SimpleMolGraph;
        vmatchgen=vmatchgen, ematchgen=ematchgen, kwargs...)
    # Note: InChI is better if you don't need mapping
    exact_match_prefilter(mol1, mol2) || return ()
    return isomorphisms(
        mol1.graph, mol2.graph,
        vmatch=vmatchgen(mol1, mol2), ematch=ematchgen(mol1, mol2); kwargs...)
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
        vmatchgen=vmatchgen, ematchgen=ematchgen, kwargs...)
    (nv(mol1) == 0 || nv(mol2) == 0) && return ()
    substruct_match_prefilter(mol1, mol2) || return ()
    return subgraph_monomorphisms(
        mol1.graph, mol2.graph,
        vmatch=vmatchgen(mol1, mol2), ematch=ematchgen(mol1, mol2); kwargs...)
end


"""
    has_substruct_match(mol1, mol2; kwargs...) -> Bool

Return whether `mol` has `query` as a substructure.
See [`substruct_matches`](@ref) for available options.
"""
has_substruct_match(mol1, mol2; kwargs...) = !isempty(substruct_matches(mol1, mol2; kwargs...))


"""
    node_substruct_matches(mol1, mol2; kwargs...) -> Iterator

Return a lazy iterator that generate node mappings between `mol` and `query` if `mol` has `query` as a substructure.
See [`substruct_matches`](@ref) for available options.
"""
function node_substruct_matches(mol1::SimpleMolGraph, mol2::SimpleMolGraph;
        vmatchgen=vmatchgen, ematchgen=ematchgen, kwargs...)
    (nv(mol1) == 0 || nv(mol2) == 0) && return ()
    substruct_match_prefilter(mol1, mol2) || return ()
    return nodesubgraph_isomorphisms(
        mol1.graph, mol2.graph,
        vmatch=vmatchgen(mol1, mol2), ematch=ematchgen(mol1, mol2); kwargs...)
end


"""
    has_node_substruct_match(mol1, mol2; kwargs...) -> Bool

Return whether `mol` has `query` as a substructure.
See [`substruct_matches`](@ref) for available options.
"""
has_node_substruct_match(mol1, mol2; kwargs...) = !isempty(node_substruct_matches(mol1, mol2; kwargs...))


"""
    edge_substruct_matches(mol1, mol2; kwargs...) -> Iterator

Return a lazy iterator that generate node mappings between `mol` and `query` if `mol` has `query` as a substructure.
See [`substruct_matches`](@ref) for available options.
"""
function edge_substruct_matches(mol1::SimpleMolGraph, mol2::SimpleMolGraph;
        vmatchgen=vmatchgen, ematchgen=ematchgen, kwargs...)
    (ne(mol1) == 0 || ne(mol2) == 0) && return ()
    substruct_match_prefilter(mol1, mol2) || return ()
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
has_edge_substruct_match(mol1, mol2; kwargs...) = !isempty(edge_substruct_matches(mol1, mol2; kwargs...))



# MCS


mcis_constraints(
    mol::SimpleMolGraph, vmatchvecgen=vmatchvecgen, ematchvecgen=ematchvecgen
) = mcis_constraints(
    Val{:connection}(), mol.graph, vmatchvec=vmatchvecgen(mol), ematchvec=ematchvecgen(mol)
)

mces_constraints(
    mol::SimpleMolGraph, vmatchvecgen=vmatchvecgen, ematchvecgen=ematchvecgen
) = mces_constraints(
    Val{:connection}(), mol.graph, vmatchvec=vmatchvecgen(mol), ematchvec=ematchvecgen(mol)
)

tdmcis_constraints(
    mol::SimpleMolGraph; vmatchvecgen=vmatchvecgen, ematchvecgen=ematchvecgen, kwargs...
) = mcis_constraints(
    Val{:shortest}(), mol.graph, vmatchvec=vmatchvecgen(mol), ematchvec=ematchvecgen(mol); kwargs...
)

tdmces_constraints(
    mol::SimpleMolGraph; vmatchvecgen=vmatchvecgen, ematchvecgen=ematchvecgen, kwargs...
) = mces_constraints(
    Val{:shortest}(), mol.graph, vmatchvec=vmatchvecgen(mol), ematchvec=ematchvecgen(mol); kwargs...
)


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
    mol1_ = mcis_constraints(mol1)
    mol2_ = mcis_constraints(mol2)
    return maximum_common_subgraph(mol1_, mol2_; kwargs...)
end

function disconnected_mces(mol1::SimpleMolGraph, mol2::SimpleMolGraph; kwargs...)
    mol1_ = mces_constraints(mol1)
    mol2_ = mces_constraints(mol2)
    return maximum_common_subgraph(mol1_, mol2_; kwargs...)
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
    mol1_ = mcis_constraints(mol1)
    mol2_ = mcis_constraints(mol2)
    return maximum_common_subgraph(mol1_, mol2_, connected=true; kwargs...)
end

function connected_mces(mol1::SimpleMolGraph, mol2::SimpleMolGraph; kwargs...)
    mol1_ = mces_constraints(mol1)
    mol2_ = mces_constraints(mol2)
    return maximum_common_subgraph(mol1_, mol2_, connected=true; kwargs...)
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
    mol1_ = tdmcis_constraints(mol1; kwargs...)
    mol2_ = tdmcis_constraints(mol2; kwargs...)
    return maximum_common_subgraph(mol1_, mol2_; kwargs...)
end

function tdmces(mol1::SimpleMolGraph, mol2::SimpleMolGraph; kwargs...)
    mol1_ = tdmces_constraints(mol1; kwargs...)
    mol2_ = tdmces_constraints(mol2; kwargs...)
    return maximum_common_subgraph(mol1_, mol2_; kwargs...)
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
