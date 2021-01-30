#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    structmatches,
    exactmatches, hasexactmatch,
    substructmatches, hassubstructmatch,
    nodeinducedmatches, hasnodeinducedmatch,
    edgeinducedmatches, hasedgeinducedmatch


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
- `fastsingleton::Bool`: if true, skip vf2 if the query is single node or edge that can be matched by simple iteration (dafault: true)
- `atommatcher::Function`: a function for semantic atom attribute matching (default: `MolecularGraph.atommatch`)
- `bondmatcher::Function`: a function for semantic bond attribute matching (default: `MolecularGraph.edgematch`)
- `mandatory::Dict{Int,Int}`: mandatory node mapping (or edge mapping if matchtype=:edgeinduced)
- `forbidden::Dict{Int,Int}`: forbidden node mapping (or edge mapping if matchtype=:edgeinduced)
- `timeout::Union{Int,Nothing}`: if specified, abort vf2 calculation when the time reached and return empty iterator (default: 10 seconds).

Note that if the query is disconnected (has component level expression like "CCO.O.O"), monomorphism mapping can be extremely slow.

Note that null mol and null query never match (e.g. isstructmatch(smilestomol(""), smilestomol("")) is false)
"""
function structmatches(
        mol::UndirectedGraph, query::UndirectedGraph, matchtype;
        prefilter=true, fastsingleton=true,
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

    # Skip isomorphism calculation if the query is a sigle node/edge
    if fastsingleton
        if nodecount(query) == 1
            matchtype !== :exact || nodecount(mol) == 1 || return ()
            res = Dict{Int,Int}[]
            q = pop!(nodeset(query))
            for n in nodeset(mol)
                if haskey(kwargs, :mandatory)
                    get(kwargs[:mandatory], n, -1) == q || continue
                end
                if haskey(kwargs, :forbidden)
                    get(kwargs[:forbidden], n, -1) == q && continue
                end
                afunc(n, q) && push!(res, Dict(n => q))
            end
            return res
        elseif nodecount(query) == 2 && edgecount(query) == 1
            matchtype !== :exact || edgecount(mol) == 1 || return ()
            res = Dict{Int,Int}[]
            q = pop!(edgeset(query))
            (qu, qv) = getedge(query, q)
            for e in edgeset(mol)
                if haskey(kwargs, :mandatory)
                    get(kwargs[:mandatory], e, -1) == q || continue
                end
                if haskey(kwargs, :forbidden)
                    get(kwargs[:forbidden], e, -1) == q && continue
                end
                bfunc(e, q) || continue
                (u, v) = getedge(mol, e)
                if afunc(u, qu) && afunc(v, qv)
                    if matchtype === :edgeinduced
                        push!(res, Dict(e => q))
                    else
                        push!(res, Dict(u => qu, v => qv))
                    end
                elseif afunc(u, qv) && afunc(v, qu)
                    if matchtype === :edgeinduced
                        push!(res, Dict(e => q))
                    else
                        push!(res, Dict(u => qv, v => qu))
                    end
                end
            end
            return res
        end
    end

    # Isomorphism
    if matchtype === :exact
        return isomorphisms(mol, query,
            nodematcher=afunc, edgematcher=bfunc; kwargs...)
    elseif matchtype === :substruct
        return subgraph_monomorphisms(mol, query,
            nodematcher=afunc, edgematcher=bfunc; kwargs...)
    elseif matchtype === :nodeinduced
        return nodesubgraph_isomorphisms(mol, query,
            nodematcher=afunc, edgematcher=bfunc; kwargs...)
    elseif matchtype === :edgeinduced
        return edgesubgraph_isomorphisms(mol, query,
            nodematcher=afunc, edgematcher=bfunc; kwargs...)
    end

    # Monomorphism by maximum cardinality matching
    # returns only the first match
    @assert matchtype === :substructmc
    if edgecount(query) != 0
        # Edge induced subgraph mapping
        for emap in edgesubgraph_isomorphisms(
                mol, query, nodematcher=afunc, edgematcher=bfunc; kwargs...)
            # Isolated node mapping
            msub = edgesubgraph(mol, Set(keys(emap)))
            qsub = edgesubgraph(query, Set(values(emap)))
            miso = setdiff(nodeset(mol), nodeset(msub))
            qiso = setdiff(nodeset(query), nodeset(qsub))
            mrem = nodesubgraph(mol, miso)
            qrem = nodesubgraph(query, qiso)
            nmap = maxcardmap(nodeset(mrem), nodeset(qrem), afunc)
            # TODO: connectivity filter
            if length(nmap) == nodecount(qrem)
                return [nmap]
            end
        end
    else
        # Isolated nodes only
        nmap = maxcardmap(nodeset(mol), nodeset(query), afunc)
        # TODO: connectivity filter
        if length(nmap) == nodecount(query)
            return [nmap]
        end
    end
    return ()
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


function atommatch(mol1::UndirectedGraph, mol2::UndirectedGraph)
    sym1 = atomsymbol(mol1)
    sym2 = atomsymbol(mol2)
    pi1 = pielectron(mol1)
    pi2 = pielectron(mol2)
    return (a1, a2) -> sym1[a1] == sym2[a2] && pi1[a1] == pi2[a2]
end

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


function bondmatch(mol1::UndirectedGraph, mol2::UndirectedGraph)
    # TODO: need bond attribute matching?
    return (b1, b2) -> true
end

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


function querymatchtree(
        query::Pair, mol::UndirectedGraph, matcher::Dict, i::Int)
    if query.first == :any
        return true
    elseif query.first == :and
        return all(querymatchtree(q, mol, matcher, i) for q in query.second)
    elseif query.first == :or
        return any(querymatchtree(q, mol, matcher, i) for q in query.second)
    elseif query.first == :not
        return !querymatchtree(query.second, mol, matcher, i)
    elseif query.first == :stereo
        # TODO: stereo not implemented yet
        return true
    elseif query.first == :recursive
        subq = parse(SMARTS, query.second)
        return hassubstructmatch(mol, subq, mandatory=Dict(i => 1))
    else
        if query.first == :sssrsizes
            return query.second in matcher[query.first][i]
        else
            return matcher[query.first][i] == query.second
        end
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
