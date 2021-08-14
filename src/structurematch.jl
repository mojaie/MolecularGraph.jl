#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    atommatch, bondmatch,
    isequivalent, query_contains, isaatommatch, isabondmatch,
    structmatches,
    exactmatches, hasexactmatch,
    substructmatches, hassubstructmatch,
    nodeinducedmatches, hasnodeinducedmatch,
    edgeinducedmatches, hasedgeinducedmatch,
    maxcommonsubstruct,
    connectedmcis, connectedmces,
    disconnectedmcis, disconnectedmces,
    tcmcis, tcmces


"""
    atommatch(mol1::UndirectedGraph, mol2::UndirectedGraph) -> Function

Return a default atom attribute comparator between two atoms.
"""
function atommatch(mol1::UndirectedGraph, mol2::UndirectedGraph)
    sym1 = atomsymbol(mol1)
    sym2 = atomsymbol(mol2)
    pi1 = pielectron(mol1)
    pi2 = pielectron(mol2)
    return (a1, a2) -> sym1[a1] == sym2[a2] && pi1[a1] == pi2[a2]
end


"""
    bondmatch(mol1::UndirectedGraph, mol2::UndirectedGraph) -> Function

Return a default bond attribute comparator between two bonds.
"""
function bondmatch(mol1::UndirectedGraph, mol2::UndirectedGraph)
    # TODO: need bond attribute matching?
    return (b1, b2) -> true
end


function querymatchtree(
        query::QueryFormula, mol::UndirectedGraph, matcher::Dict, i::Int)
    if query.key === :any
        return true
    elseif query.key === :and
        return all(querymatchtree(q, mol, matcher, i) for q in query.value)
    elseif query.key === :or
        return any(querymatchtree(q, mol, matcher, i) for q in query.value)
    elseif query.key === :not
        return !querymatchtree(query.value, mol, matcher, i)
    elseif query.key === :stereo
        # TODO: stereo not implemented yet
        return true
    elseif query.key === :recursive
        subq = parse(SMARTS, query.value)
        return hassubstructmatch(mol, subq, mandatory=Dict(i => 1))
    else
        if query.key === :sssrsizes
            return query.value in matcher[query.key][i]
        else
            return matcher[query.key][i] == query.value
        end
    end
end


"""
    atommatch(mol::UndirectedGraph, querymol::QueryMol) -> Function

Return a default atom attribute comparator that returns true if the atom satisfies all the queryatom conditions.
"""
function atommatch(mol::UndirectedGraph, querymol::QueryMol)
    matcher = Dict(
        :atomsymbol => atomsymbol(mol),
        :isaromatic => isaromatic(mol),
        :charge => charge(mol),
        :mass => getproperty.(nodeattrs(mol), :mass),
        :stereo => getproperty.(nodeattrs(mol), :stereo),
        :connectivity => connectivity(mol),
        :nodedegree => nodedegree(mol),
        :valence => valence(mol),
        :hydrogenconnected => hydrogenconnected(mol),
        :sssrsizes => sssrsizes(mol),
        :sssrcount => sssrcount(mol)
    )
    return (a, qa) -> querymatchtree(
        nodeattr(querymol, qa).query, mol, matcher, a)
end


"""
    bondmatch(mol::UndirectedGraph, querymol::QueryMol) -> Function

Return a default bond attribute comparator that returns true if the bond satisfies all the querybond conditions.
"""
function bondmatch(mol::UndirectedGraph, querymol::QueryMol)
    matcher = Dict(
        :bondorder => bondorder(mol),
        :isringbond => isringbond(mol),
        :isaromaticbond => isaromaticbond(mol),
        :stereo => getproperty.(edgeattrs(mol), :stereo)
    )
    return (b, qb) -> querymatchtree(
        edgeattr(querymol, qb).query, mol, matcher, b)
end


"""
    exactatommatch(qmol1::QueryMol, qmol2::QueryMol) -> Function

Return a default atom attribute comparator between two atom queries.
"""
function exactatommatch(qmol1::QueryMol, qmol2::QueryMol)
    return function (a1, a2)
        return nodeattr(qmol1, a1).query == nodeattr(qmol2, a2).query
    end
end


"""
    exactbondmatch(qmol1::QueryMol, qmol2::QueryMol) -> Function

Return a default bond attribute comparator between two bond queries.
"""
function exactbondmatch(qmol1::QueryMol, qmol2::QueryMol)
    return function (b1, b2)
        return edgeattr(qmol1, b1).query == edgeattr(qmol2, b2).query
    end
end


"""
    atommatch(qmol1::QueryMol, qmol2::QueryMol) -> Function

Return an atom attribute comparator that returns true if a2 contains a1.
"""
function atommatch(qmol1::QueryMol, qmol2::QueryMol)
    return function (a1, a2)
        return issubset(nodeattr(qmol1, a1).query, nodeattr(qmol2, a2).query)
    end
end


"""
    bondmatch(qmol1::QueryMol, qmol2::QueryMol) -> Function

Return a bond attribute comparator that returns true if b2 contains b1.
"""
function bondmatch(qmol1::QueryMol, qmol2::QueryMol)
    return function (b1, b2)
        return issubset(edgeattr(qmol1, b1).query, edgeattr(qmol2, b2).query)
    end
end



function fastidentityfilter(mol1::UndirectedGraph, mol2::UndirectedGraph)
    if nodecount(mol1) != nodecount(mol2)
        return false
    elseif edgecount(mol1) != edgecount(mol2)
        return false
    elseif circuitrank(mol1) != circuitrank(mol2)
        return false
    end
    return true
end


function fastsubstrfilter(mol1::UndirectedGraph, mol2::UndirectedGraph)
    if nodecount(mol1) < nodecount(mol2)
        return false
    elseif edgecount(mol1) < edgecount(mol2)
        return false
    elseif circuitrank(mol1) < circuitrank(mol2)
        return false
    end
    return true
end
    

"""
    structmatches(
        mol::UndirectedGraph, query::UndirectedGraph,
        matchtype::Symbol; kwargs...) -> Iterator

Return a lazy iterator that generate all isomorphism mappings between `mol` and `query`.

`matchtype` should be one of the followings
  - `:exact`: mol is isomorphic to query
  - `:substruct`: a subgraph of mol is monomorphic to query (typical substructure search)
  - `:nodeinduced`: a node-induced subgraph of mol is isomorphic to query
  - `:edgeinduced`: an edge-induced subgraph of mol is isomorphic to query`

# options

- `prefilter::Bool`: if true, apply simple prefilter by graph size and topology to skip vf2 calculation (dafault: true)
- `fastsingleton::Bool`: if true, skip vf2 if the query is single node or edge that can be matched by simple iteration (dafault: false)
- `atommatcher::Function`: a function for semantic atom attribute matching (default: `MolecularGraph.atommatch`)
- `bondmatcher::Function`: a function for semantic bond attribute matching (default: `MolecularGraph.edgematch`)
- `mandatory::Dict{Int,Int}`: mandatory node mapping (or edge mapping if matchtype=:edgeinduced)
- `forbidden::Dict{Int,Int}`: forbidden node mapping (or edge mapping if matchtype=:edgeinduced)
- `timeout::Union{Int,Nothing}`: if specified, abort vf2 calculation when the time reached and return empty iterator (default: 10 seconds).

Note that null mol and null query never match (e.g. isstructmatch(smilestomol(""), smilestomol("")) is false)
"""
function structmatches(
        mol::UndirectedGraph, query::UndirectedGraph, matchtype;
        prefilter=true, fastsingleton=false,
        atommatcher=atommatch, bondmatcher=bondmatch, kwargs...)
    # Null molecule filter
    nodecount(mol) == 0 && return ()
    nodecount(query) == 0 && return ()
    matchtype === :edgeinduced && edgecount(mol) == 0 && return ()
    matchtype === :edgeinduced && edgecount(query) == 0 && return ()

    # Prefilter by graph size and topology
    if prefilter
        if matchtype === :exact
            fastidentityfilter(mol, query) || return ()
        else
            fastsubstrfilter(mol, query) || return ()
        end
    end

    # Matcher functions
    afunc = atommatcher(mol, query)
    bfunc = bondmatcher(mol, query)

    # Skip substruct calculation if the query is a sigle node/edge
    # This slightly speeds up SMARTS queries and hence functional group analysis
    # Note that mandatory and forbidden mapping would be ignored
    if matchtype === :substruct && fastsingleton
        if nodecount(query) == 1
            res = Dict{Int,Int}[]
            q = pop!(nodeset(query))
            for n in nodeset(mol)
                afunc(n, q) && push!(res, Dict(n => q))
            end
            return res
        elseif nodecount(query) == 2 && edgecount(query) == 1
            res = Dict{Int,Int}[]
            q = pop!(edgeset(query))
            (qu, qv) = getedge(query, q)
            for e in edgeset(mol)
                bfunc(e, q) || continue
                (u, v) = getedge(mol, e)
                if afunc(u, qu) && afunc(v, qv)
                    push!(res, Dict(u => qu, v => qv))
                elseif afunc(u, qv) && afunc(v, qu)
                    push!(res, Dict(u => qv, v => qu))
                end
            end
            return res
        end
    end

    # Isomorphism
    if matchtype === :exact
        return isomorphisms(mol, query, nodematcher=afunc, edgematcher=bfunc; kwargs...)
    elseif matchtype === :substruct
        return subgraph_monomorphisms(mol, query, nodematcher=afunc, edgematcher=bfunc; kwargs...)
    elseif matchtype === :nodeinduced
        return nodesubgraph_isomorphisms(mol, query, nodematcher=afunc, edgematcher=bfunc; kwargs...)
    elseif matchtype === :edgeinduced
        return edgesubgraph_isomorphisms(mol, query, nodematcher=afunc, edgematcher=bfunc; kwargs...)
    end
end


"""
    exactmatches(
        mol::UndirectedGraph, query::UndirectedGraph; kwargs...) -> Iterator

Return a lazy iterator that generate node mappings between `mol` and `query` if these are exactly same.
See [`MolecularGraph.structmatches`](@ref) for available options.
"""
exactmatches(mol1, mol2; kwargs...
    ) = structmatches(mol1, mol2, :exact; kwargs...)


"""
    hasexactmatch(
        mol::UndirectedGraph, query::UndirectedGraph; kwargs...) -> Iterator

Return whether `mol` and `query` have exactly the same structure.
See [`MolecularGraph.structmatches`](@ref) for available options.
"""
hasexactmatch(mol1, mol2; kwargs...
    ) = !isempty(exactmatches(mol1, mol2; kwargs...))


"""
    substructmatches(
        mol::UndirectedGraph, query::UndirectedGraph; kwargs...) -> Iterator

Return a lazy iterator that generate node mappings between `mol` and `query` if `mol` has `query` as a substructure.
See [`MolecularGraph.structmatches`](@ref) for available options.
"""
substructmatches(mol1, mol2; kwargs...
    ) = structmatches(mol1, mol2, :substruct; kwargs...)


"""
    hassubstructmatch(
        mol::UndirectedGraph, query::UndirectedGraph; kwargs...) -> Iterator

Return whether `mol` has `query` as a substructure.
See [`MolecularGraph.structmatches`](@ref) for available options.
"""
hassubstructmatch(mol1, mol2; kwargs...
    ) = !isempty(substructmatches(mol1, mol2; kwargs...))


"""
    nodeinducedmatches(
        mol::UndirectedGraph, query::UndirectedGraph; kwargs...) -> Iterator

Return a lazy iterator that generate node mappings between `mol` and `query` if a node-induced subgraph of `mol` graph is isomorphic to `query` graph.
See [`MolecularGraph.structmatches`](@ref) for available options.
"""
nodeinducedmatches(mol1, mol2; kwargs...
    ) = structmatches(mol1, mol2, :nodeinduced; kwargs...)


"""
    hasnodeinducedmatch(
        mol::UndirectedGraph, query::UndirectedGraph; kwargs...) -> Iterator

Return whether a node-induced subgraph of `mol` graph is isomorphic to `query` graph.
See [`MolecularGraph.structmatches`](@ref) for available options.
"""
hasnodeinducedmatch(mol1, mol2; kwargs...
    ) = !isempty(nodeinducedmatches(mol1, mol2; kwargs...))


"""
    edgeinducedmatches(
        mol::UndirectedGraph, query::UndirectedGraph; kwargs...) -> Iterator

Return a lazy iterator that generate **edge** mappings between `mol` and `query` if a edge-induced subgraph of `mol` graph is isomorphic to `query` graph.
See [`MolecularGraph.structmatches`](@ref) for available options.
"""
edgeinducedmatches(mol1, mol2; kwargs...
    ) = structmatches(mol1, mol2, :edgeinduced; kwargs...)


"""
    hasedgeinducedmatch(
        mol::UndirectedGraph, query::UndirectedGraph; kwargs...) -> Iterator

Return whether a edge-induced subgraph of `mol` graph is isomorphic to `query` graph.
See [`MolecularGraph.structmatches`](@ref) for available options.
"""
hasedgeinducedmatch(mol1, mol2; kwargs...
    ) = !isempty(edgeinducedmatches(mol1, mol2; kwargs...))



"""
    maxcommonsubstruct(mol1::UndirectedGraph, mol2::UndirectedGraph, matchtype; kwargs...) -> MaxCommonSubgraphResult

Compute maximum common substructure (MCS) of mol1 and mol2.

## Keyword arguments

- matchtype(Symbol): :nodeinduced or :edgeinduced (default=:edgeinduced).
- connected(Bool): if true, apply connected MCS constraint (default=false).
- topological(Bool): if true, apply topological constraint (default=false).
- diameter(Int): distance cutoff for topological constraint.
- tolerance(Int): distance mismatch tolerance for topological constraint.
- timeout(Int): abort calculation and return suboptimal results if the execution
time has reached the given value (default=60, in seconds).
- targetsize(Int): abort calculation and return suboptimal result so far if the
given mcs size achieved.

# References

1. Kawabata, T. (2011). Build-Up Algorithm for Atomic Correspondence between
Chemical Structures. Journal of Chemical Information and Modeling, 51(8),
1775–1787. https://doi.org/10.1021/ci2001023
1. https://www.jstage.jst.go.jp/article/ciqs/2017/0/2017_P4/_article/-char/en
"""
function maxcommonsubstruct(mol1::UndirectedGraph, mol2::UndirectedGraph, matchtype; kwargs...)
    afunc = atommatch(mol1, mol2)
    bfunc = bondmatch(mol1, mol2)
    return maxcommonsubgraph(mol1, mol2, matchtype, nodematcher=afunc, edgematcher=bfunc; kwargs...)
end

disconnectedmcis(mol1, mol2; kwargs...) = maxcommonsubstruct(mol1, mol2, :nodeinduced; kwargs...)

disconnectedmces(mol1, mol2; kwargs...) = maxcommonsubstruct(mol1, mol2, :edgeinduced; kwargs...)

connectedmcis(mol1, mol2; kwargs...) = maxcommonsubstruct(mol1, mol2, :nodeinduced, connected=true; kwargs...)

connectedmces(mol1, mol2; kwargs...) = maxcommonsubstruct(mol1, mol2, :edgeinduced, connected=true; kwargs...)

tcmcis(mol1, mol2; kwargs...) = maxcommonsubstruct(mol1, mol2, :nodeinduced, topological=true; kwargs...)

tcmces(mol1, mol2; kwargs...) = maxcommonsubstruct(mol1, mol2, :edgeinduced, topological=true; kwargs...)



"""
    nmap = emaptonmap(emap, mol, query)

Convert an edge-based mapping, of the form returned by [`edgesubgraphmatches`](@ref), into
a map between nodes. Commonly, `nmap[i]` is a length-1 vector `[j]`, where `i=>j` is the mapping
from `nodeattr(query, i)` to `nodeattr(mol, j)`. In cases where the mapping is ambiguous,
`nmap[i]` may be multivalued.
"""
function emaptonmap(emap, mol::UndirectedGraph, query::QueryMol)
    nmol, nq = nodecount(mol), nodecount(query)
    nq <= nmol || throw(ArgumentError("query must be a substructure of mol"))
    # Each node in the query edges can map to either of two nodes in mol
    qnodeoptions = [Tuple{Int,Int}[] for _ = 1:nq]
    for (edgeidmol, edgeidq) in emap
        edgemol, edgeq = getedge(mol, edgeidmol), getedge(query, edgeidq)
        for nq in edgeq
            push!(qnodeoptions[nq], edgemol)
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
