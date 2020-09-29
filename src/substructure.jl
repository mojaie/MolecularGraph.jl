#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    structmatches, isstructmatch, issmartsgroupmatch


"""

    structmatches(mol::UndirectedGraph, query::UndirectedGraph; kwargs...) -> Iterator

query is a subgraph of mol


Generate molecular graph match mappings between `mol1` and `mol2`. If no match
found, return nothing.

option
matchtype: exact, substruct, nodeinduced, edgeinduced

prefilter(true): if true, apply prefilter before vf2
fastsingleton(true): skip vf2 if the query is single node or edge

nodematcher(atommatch(mol1,mol2))
edgematcher(edgematch(mol1,mol2))
mandatory(nothing): invalid if induced=none
forbidden(nothihg): invalid if induced=none
timeout(10): if specified, set timeout

smarts: induced=none

The mapping is based on only edge induced subgraph isomorphism and therefore
it does not care disconnected single atom matches. This function is intended
for use in substructure search. If you need single atom SMARTS match
(ex. [#16;X2;!R]), see [`querymatch`](@ref).


    structmatches(mol1::UndirectedGraph, mol2::UndirectedGraph; kwargs...) -> Iterator

Generate substructure match mappings between `mol1` and `mol2`. If no match
found, return nothing.



querymatch(mol::UndirectedGraph, query::QueryMol; kwargs...) -> Dict{Int,Int}

Generate substructure match mappings between `mol1` and `mol2`. If no match
found, return nothing.

This accepts also disconnected single atom but returns only the first match.
This function is intended for use in SMARTS query search

fastquerymatches(mol::UndirectedGraph, query::QueryMol; kwargs...
    ) -> Dict{Int,Int}

Generate query match mappings between `mol` and `query`. If no match
found, return nothing.

The `query` should not have any component level expression that means it should
not have any dots (`.`). This is intended for use in functional group detection.

Note that null mol and null query never match
(e.g. isstructmatch(smilestomol(""), smilestomol("")) is false)
"""
function structmatches(
        mol::UndirectedGraph, query::UndirectedGraph, matchtype;
        prefilter=true, fastsingleton=true,
        atommatcher=atommatch, bondmatcher=bondmatch, kwargs...)
    # Null molecule
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
        return exactmatches(mol, query,
            nodematcher=afunc, edgematcher=bfunc; kwargs...)
    elseif matchtype === :nodeinduced
        return subgraphmatches(mol, query,
            nodematcher=afunc, edgematcher=bfunc; kwargs...)
    elseif matchtype === :edgeinduced
        return edgesubgraphmatches(mol, query,
            nodematcher=afunc, edgematcher=bfunc; kwargs...)
    end

    # TODO: Monomorphism by maximum cardinality mapping
    # returns only the first match
    @assert matchtype === :substruct
    if edgecount(query) != 0
        # Edge induced subgraph mapping
        for emap in edgesubgraphmatches(
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


isstructmatch(mol1, mol2, matchtype; kwargs...
    ) = !isempty(structmatches(mol1, mol2, matchtype; kwargs...))



function issmartsgroupmatch(
        mol::UndirectedGraph, query::QueryMol, root::Int; kwargs...)
    # For recursive SMARTS match
    isempty(query.connectivity) || throw(
        ErrorException("Component level query is disallowed"))
    afunc = atommatch(mol, query)
    bfunc = bondmatch(mol, query)
    @assert root in nodeset(mol) "node $(root) does not exist"
    if nodecount(query) == 1
        # node match
        return afunc(root, 1)
    elseif edgecount(query) == 1
        # edge match
        for (inc, adj) in neighbors(mol, root)
            if afunc(root, 1) && afunc(adj, 2) && bfunc(inc, 1)
                return true
            end
        end
        return false
    else
        # subgraph match
        for n in incidences(mol, root)
            if isedgesubgraphmatch(
                    mol, query, nodematcher=afunc, edgematcher=bfunc,
                    mandatory=Dict(n=>1); kwargs...)
                return true
            end
        end
        return false
    end
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
        return issmartsgroupmatch(mol, subq, i)
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
