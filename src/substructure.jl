#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    structmatches, structmatch, isstructmatch,
    substructmatches, substructmatch, issubstructmatch,
    querymatch, isquerymatch,
    fastquerymatches, isqueryidentical, issmartsgroupmatch


"""
    structmatches(mol1::UndirectedGraph, mol2::UndirectedGraph; kwargs...) -> Iterator

Generate molecular graph match mappings between `mol1` and `mol2`. If no match
found, return nothing.
"""
function structmatches(mol1::UndirectedGraph, mol2::UndirectedGraph;
        prefilter=true, kwargs...)
    !prefilter || fastidentityfilter(mol1, mol2) || return []
    afunc = atommatch(mol1, mol2)
    bfunc = bondmatch(mol1, mol2)
    return graphmatches(mol1, mol2,
        nodematcher=afunc, edgematcher=bfunc; kwargs...)
end

function structmatch(mol1, mol2; kwargs...)
    res = iterate(structmatches(mol1, mol2; kwargs...))
    return res === nothing ? nothing : res[1]
end

isstructmatch(mol1, mol2; kwargs...
    ) = structmatch(mol1, mol2; kwargs...) !== nothing


"""
    substructmatches(mol1::UndirectedGraph, mol2::UndirectedGraph; kwargs...) -> Iterator

Generate substructure match mappings between `mol1` and `mol2`. If no match
found, return nothing.

The mapping is based on only edge induced subgraph isomorphism and therefore
it does not care disconnected single atom matches. This function is intended
for use in substructure search. If you need single atom SMARTS match
(ex. [#16;X2;!R]), see [`querymatch`](@ref).
"""
function substructmatches(mol1::UndirectedGraph, mol2::UndirectedGraph;
        prefilter=true, kwargs...)
    !prefilter || fastsubstrfilter(mol1, mol2) || return []
    afunc = atommatch(mol1, mol2)
    bfunc = bondmatch(mol1, mol2)
    return edgesubgraphmatches(mol1, mol2,
        nodematcher=afunc, edgematcher=bfunc; kwargs...)
end

function substructmatch(mol1, mol2; kwargs...)
    res = iterate(substructmatches(mol1, mol2; kwargs...))
    return res === nothing ? nothing : res[1]
end

issubstructmatch(mol1, mol2; kwargs...
    ) = substructmatch(mol1, mol2; kwargs...) !== nothing



"""
    querymatch(mol::UndirectedGraph, query::QueryMol; kwargs...) -> Dict{Int,Int}

Generate substructure match mappings between `mol1` and `mol2`. If no match
found, return nothing.

This accepts also disconnected single atom but returns only the first match.
This function is intended for use in SMARTS query search

"""
function querymatch(mol::UndirectedGraph, query::QueryMol; kwargs...)
    mnodeset = nodeset(mol)
    qnodeset = nodeset(query)
    afunc = atommatch(mol, query)
    bfunc = bondmatch(mol, query)
    if edgecount(query) != 0
        # Edge induced subgraph mapping
        for emap in edgesubgraphmatches(
                mol, query, nodematcher=afunc, edgematcher=bfunc; kwargs...)
            # Isolated node mapping
            msub = edgesubgraph(mol, Set(keys(emap)))
            qsub = edgesubgraph(query, Set(values(emap)))
            miso = setdiff(mnodeset, nodeset(msub))
            qiso = setdiff(qnodeset, nodeset(qsub))
            mrem = nodesubgraph(mol, miso)
            qrem = nodesubgraph(query, qiso)
            nmap = maxcardmap(nodeset(mrem), nodeset(qrem), afunc)
            # TODO: connectivity filter
            if length(nmap) == nodecount(qrem)
                return (emap, nmap)
            end
        end
    elseif nodecount(query) != 0
        # Isolated nodes only
        nmap = maxcardmap(mnodeset, qnodeset, afunc)
        # TODO: connectivity filter
        if length(nmap) == nodecount(query)
            return (Dict(), nmap)
        end
    end
    return
end

"""
    isquerymatch(mol, query; kwargs...)

Return whether mol matches with the query.
"""
isquerymatch(mol, query; kwargs...
    ) = querymatch(mol, query; kwargs...) !== nothing


"""
    fastquerymatches(mol::UndirectedGraph, query::QueryMol; kwargs...
        ) -> Dict{Int,Int}

Generate query match mappings between `mol` and `query`. If no match
found, return nothing.

The `query` should not have any component level expression that means it should
not have any dots (`.`). This is intended for use in functional group detection.

"""
function fastquerymatches(mol::UndirectedGraph, query::QueryMol; kwargs...)
    isempty(query.connectivity) || throw(
        ErrorException("Component level query is disallowed"))
    afunc = atommatch(mol, query)
    bfunc = bondmatch(mol, query)
    nodecount(query) == 0 && return ()
    if nodecount(query) == 1
        # node match
        match = Iterators.filter(nodeset(mol)) do ma
            return afunc(ma, 1)
        end
        return ((Dict(), Dict(ma => 1)) for ma in match)
    elseif edgecount(query) == 1
        # edge match
        (qu, qv) = getedge(query, 1)
        match = Iterators.filter(edgeset(mol)) do mb
            (u, v) = getedge(mol, mb)
            return (
                ((afunc(u, qu) && afunc(v, qv))
                || (afunc(u, qv) && afunc(v, qu)))
                && bfunc(mb, 1)
            )
        end
        return ((Dict(mb => 1), Dict()) for mb in match)
    else
        # subgraph match
        return ((emap, Dict()) for emap in edgesubgraphmatches(
            mol, query, nodematcher=afunc, edgematcher=bfunc; kwargs...))
    end
end


function isqueryidentical(mol::UndirectedGraph, query::QueryMol; kwargs...)
    isempty(query.connectivity) || throw(
        ErrorException("Component level query is disallowed"))
    afunc = atommatch(mol, query)
    bfunc = bondmatch(mol, query)
    nodecount(query) == 0 && return ()
    if nodecount(query) == 1
        for n in nodeset(mol)
            afunc(n, 1) && return true
        end
        return false
    elseif edgecount(query) == 1
        (qu, qv) = getedge(query, 1)
        for e in edgeset(mol)
            (u, v) = getedge(mol, e)
            pred = ((afunc(u, qu) && afunc(v, qv))
                || (afunc(u, qv) && afunc(v, qu))) && bfunc(e, 1)
            pred && return true
        end
        return false
    else
        return isgraphmatch(
            mol, query, nodematcher=afunc, edgematcher=bfunc; kwargs...)
    end
end


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
