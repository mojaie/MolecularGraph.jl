#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    is_identical,
    is_substruct,
    is_superstruct,
    isomorphism,
    is_querymatch,
    querymatch,
    querymatchiter,
    isSMARTSgroupmatch,
    atommatch,
    bondmatch,
    fast_identity_filter,
    fast_substr_filter,
    deltaYmap!


"""
    is_identical(mol1::VectorMol, mol2::VectorMol)

Return whether mol1 and mol2 are identical in chemical structure.
"""
function is_identical(mol1::VectorMol, mol2::VectorMol)
    if !fast_identity_filter(mol1, mol2)
        return false
    end
    return isomorphism(mol1, mol2, mode=:Isomorphism) !== nothing
end


"""
    is_substruct(mol1::VectorMol, mol2::VectorMol)

Return whether mol1 is a substructure of mol2.
"""
function is_substruct(mol1::VectorMol, mol2::VectorMol)
    if !fast_substr_filter(mol2, mol1)
        return false
    end
    return isomorphism(mol2, mol1) !== nothing
end


"""
    is_superstruct(mol1, mol2)

Return whether mol1 is a superstructure of mol2.
"""
is_superstruct(mol1, mol2) = is_substruct(mol2, mol1)


function isomorphism(mol1::VectorMol, mol2::VectorMol;
        mode=:Subgraph, atommatcher=atommatch, bondmatcher=bondmatch, kwargs...)
    # Mapping based on edge subgraph isomorphism (ignore disconnected atom)
    # Intended for use in substructure search
    return edgeisomorphism(mol1.graph, mol2.graph, mode=mode,
        nodematcher=atommatcher(mol1, mol2),
        edgematcher=bondmatcher(mol1, mol2), kwargs...)
end


"""
    is_querymatch(mol, query; kwargs...)

Return whether mol matches with the query.
"""
function is_querymatch(mol, query; kwargs...)
    return querymatch(mol, query; kwargs...) !== nothing
end


function querymatch(mol::MolGraph, query::QueryMolGraph;
        atommatcher=atommatch, bondmatcher=bondmatch, kwargs...)
    # Accept also disconnected atom but return only the first match
    # Intended for use in SMARTS query search
    mnodeset = nodekeys(mol.graph)
    qnodeset = nodekeys(query.graph)
    afunc = atommatcher(mol, query)
    bfunc = bondmatcher(mol, query)
    if bondcount(query) != 0
        # Edge induced subgraph mapping
        for emap in edgeisomorphismiter(
                mol.graph, query.graph, mode=:Subgraph,
                nodematcher=afunc, edgematcher=bfunc; kwargs...)
            # Isolated node mapping
            msub = edgesubgraph(mol.graph, keys(emap))
            qsub = edgesubgraph(query.graph, values(emap))
            miso = setdiff(mnodeset, nodekeys(msub))
            qiso = setdiff(qnodeset, nodekeys(qsub))
            mrem = nodesubgraph(mol.graph, miso)
            qrem = nodesubgraph(query.graph, qiso)
            nmap = maxcardmap(nodekeys(mrem), nodekeys(qrem), afunc)
            if length(nmap) == nodecount(qrem)
                return (emap, nmap)
            end
        end
    elseif atomcount(query) != 0
        # Isolated nodes only
        nmap = maxcardmap(mnodeset, qnodeset, afunc)
        if length(nmap) == atomcount(query)
            return (Dict(), nmap)
        end
    end
    return
end


querymatch(
    mol::MolGraph, query::ConnectedQueryMol; kwargs...
) = iterate(querymatchiter(mol, query; kwargs...))


function querymatchiter(mol::MolGraph, query::ConnectedQueryMol;
        atommatcher=atommatch, bondmatcher=bondmatch, kwargs...)
    # Iterate over all possible subgraph isomorphism mappings
    # Intended for use in functional group annotation
    afunc = atommatcher(mol, query)
    bfunc = bondmatcher(mol, query)
    if atomcount(query) == 0
        return ()
    elseif atomcount(query) == 1
        # node match
        qa = pop!(nodekeys(query.graph))
        match = Iterators.filter(nodekeys(mol.graph)) do ma
            return afunc(ma, qa)
        end
        return ((Dict(), Dict(ma => qa)) for ma in match)
    elseif bondcount(query) == 1
        # edge match
        qb = pop!(edgekeys(query.graph))
        qbond = getbond(query, qb)
        match = Iterators.filter(edgekeys(mol.graph)) do mb
            mbond = getbond(mol, mb)
            return (
                ((afunc(mbond.u, qbond.u) && afunc(mbond.v, qbond.v))
                || (afunc(mbond.u, qbond.v) && afunc(mbond.v, qbond.u)))
                && bfunc(mb, qb)
            )
        end
        return ((Dict(mb => qb), Dict()) for mb in match)
    else
        # subgraph match
        return ((emap, Dict()) for emap in edgeisomorphismiter(
            mol.graph, query.graph, mode=:Subgraph,
            nodematcher=afunc, edgematcher=bfunc; kwargs...
        ))
    end
end


function isSMARTSgroupmatch(mol::MolGraph, query::ConnectedQueryMol, root;
                            atommatcher=atommatch, bondmatcher=bondmatch,
                            kwargs...)
    # For recursive SMARTS match
    afunc = atommatcher(mol, query)
    bfunc = bondmatcher(mol, query)
    @assert root in nodekeys(mol.graph) "node $(root) does not exist"
    if atomcount(query) == 1
        # node match
        return afunc(root, 1)
    elseif bondcount(query) == 1
        # edge match
        for (nbr, b) in neighbors(mol.graph, root)
            if afunc(root, 1) && afunc(nbr, 2) && bfunc(b, 1)
                return true
            end
        end
        return false
    else
        # subgraph match
        for n in neighboredgekeys(mol.graph, root)
            if is_edge_subgraph(
                    query.graph, mol.graph,
                    nodematcher=afunc, edgematcher=bfunc,
                    mandatory=Dict(n=>1); kwargs...
                )
                return true
            end
        end
        return false
    end
end


function atommatch(mol1::VectorMol, mol2::VectorMol)
    return function (a1, a2)
        sym1 = mol1[:Symbol]
        sym2 = mol2[:Symbol]
        pi1 = mol1[:Pi]
        pi2 = mol2[:Pi]
        sym1[a1] == sym2[a2] && pi1[a1] == pi2[a2]
    end
end

function atommatch(mol::VectorMol, querymol::QueryMolGraph)
    return function (a, qa)
        q = getatom(querymol, qa).query
        return querymatchtree(q, mol, a)
    end
end

atommatch(
    view::MolGraphView{G,M}, querymol::QueryMolGraph
) where {G<:SubgraphView,M<:VectorMol} = atommatch(view.molecule, querymol)



function bondmatch(mol1::VectorMol, mol2::VectorMol)
    # TODO: need bond attribute matching?
    return (b1, b2) -> true
end

function bondmatch(mol::VectorMol, querymol::QueryMolGraph)
    return function (b, qb)
        q = getbond(querymol, qb).query
        return querymatchtree(q, mol, b)
    end
end

bondmatch(
    view::MolGraphView{G,M}, querymol::QueryMolGraph
) where {G<:SubgraphView,M<:VectorMol} = bondmatch(view.molecule, querymol)


function querymatchtree(query::Pair, mol::VectorMol, i::Int)
    if query.first == :any
        return true
    elseif query.first == :and
        return all(querymatchtree(q, mol, i) for q in query.second)
    elseif query.first == :or
        return any(querymatchtree(q, mol, i) for q in query.second)
    elseif query.first == :not
        return !querymatchtree(query.second, mol, i)
    elseif query.first == :stereo
        # TODO: stereo not implemented yet
        return true
    elseif query.first == :recursive
        subq = parse(ConnectedSMARTS, query.second)
        return isSMARTSgroupmatch(mol, subq, i)
    else
        if query.first == :RingSize
            # TODO:
            return query.second in mol[query.first][i]
        else
            return mol[query.first][i] == query.second
        end
    end
end


function fast_identity_filter(mol1::VectorMol, mol2::VectorMol)
    if atomcount(mol1) != atomcount(mol2)
        return false
    elseif bondcount(mol1) != bondcount(mol2)
        return false
    elseif (length(mol1.annotation[:Topology].rings)
            != length(mol2.annotation[:Topology].rings))
        return false
    end
    return true
end


function fast_substr_filter(mol1::VectorMol, mol2::VectorMol)
    if atomcount(mol1) < atomcount(mol2)
        return false
    elseif bondcount(mol1) < bondcount(mol2)
        return false
    elseif (length(mol1.annotation[:Topology].rings)
            < length(mol2.annotation[:Topology].rings))
        return false
    end
    return true
end
