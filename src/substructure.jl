#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    molidentstate,
    molsubstrstate,
    substructmap!,
    is_identical,
    is_substruct,
    is_superstruct,
    match_molquery,
    match_groupquery,
    preprocess!,
    atommatch,
    bondmatch,
    fast_identity_filter,
    fast_substr_filter,
    deltaYmap!


molidentstate(G, H, nmatch, ematch) = VF2EdgeInducedState(
    G, H, :isomorphic, 1000, nmatch, ematch, Dict(), Dict(),
    Dict(), Dict(), Dict(), Dict(), []
)

molsubstrstate(G, H, nmatch, ematch) = VF2EdgeInducedState(
    G, H, :subgraph, 1000, nmatch, ematch, Dict(), Dict(),
    Dict(), Dict(), Dict(), Dict(), []
)


function substructmap!(mol, query, state)
    # TODO: CQS
    # TODO: not efficient in the case of small query (one atom or one bond)
    # TODO: VF2 lazy iterator
    substrs = []
    mnodeset = nodekeys(mol.graph)
    qnodeset = nodekeys(query.graph)
    matchfunc = atommatch(mol, query)
    if bondcount(query) != 0
        # Edge induced subgraph mapping
        isomorphmap!(state)
        for emap in state.mappings
            # Isolated node mapping
            msub = edgesubgraph(mol.graph, keys(emap))
            qsub = edgesubgraph(query.graph, values(emap))
            miso = setdiff(mnodeset, nodekeys(msub))
            qiso = setdiff(qnodeset, nodekeys(qsub))
            mrem = nodesubgraph(mol.graph, miso)
            qrem = nodesubgraph(query.graph, qiso)
            nmap = maxcardmap(nodekeys(mrem), nodekeys(qrem), matchfunc)
            if length(nmap) == nodecount(qrem)
                push!(substrs, (emap, nmap))
            end
        end
    elseif atomcount(query) != 0
        # Isolated nodes only
        nmap = maxcardmap(mnodeset, qnodeset, matchfunc)
        if length(nmap) == atomcount(query)
            push!(substrs, ([], nmap))
        end
    end
    return substrs
end


function is_identical(mol1::VectorMol, mol2::VectorMol)
    preprocess!(mol1)
    preprocess!(mol2)
    if !fast_identity_filter(mol1, mol2)
        return false
    end
    state = molidentstate(
        mol1.graph, mol2.graph, atommatch(mol1, mol2), bondmatch(mol1, mol2))
    return !isempty(substructmap!(mol1, mol2, state))
end


function is_substruct(mol1::VectorMol, mol2::VectorMol)
    preprocess!(mol1)
    preprocess!(mol2)
    if !fast_substr_filter(mol2, mol1)
        return false
    end
    state = molsubstrstate(
        mol2.graph, mol1.graph, atommatch(mol1, mol2), bondmatch(mol1, mol2))
    return !isempty(substructmap!(mol2, mol1, state))
end

is_superstruct(mol1, mol2) = is_substruct(mol2, mol1)


function match_molquery(mol::VectorMol, query::QueryMol)
    preprocess!(mol)
    state = molsubstrstate(
        mol.graph, query.graph, atommatch(mol, query), bondmatch(mol, query))
    return !isempty(substructmap!(mol, query, state))
end


function match_groupquery(mol::VectorMol, query::QueryMol, root::Int)
    preprocess!(mol)
    for nbr in neighboredgekeys(mol.graph, root)
        state = molsubstrstate(
            mol.graph, query.graph,
            atommatch(mol, query), bondmatch(mol, query))
        state.mandatory[nbr] = 1
        isomorphmap!(state)
        if !isempty(state.mappings)
            return true
        end
    end
    return false
end


function preprocess!(mol)
    # TODO: remove and annotate Hs and trivials
    required_annotation(mol, :Topology)
    required_annotation(mol, :Elemental)
    return
end


function atommatch(mol1::VectorMol, mol2::VectorMol)
    return function (a1, a2)
        sym1 = mol1.v[:Symbol]
        sym2 = mol2.v[:Symbol]
        pi1 = mol1.v[:Pi]
        pi2 = mol2.v[:Pi]
        sym1[a1] == sym2[a2] && pi1[a1] == pi2[a2]
    end
end

function atommatch(mol::VectorMol, querymol::QueryMol)
    return function (a, qa)
        q = getatom(querymol, qa).query
        return querymatchtree(q, mol, a)
    end
end

atommatch(query::QueryMol, mol::VectorMol) = atom_match(mol, query)


function bondmatch(mol1::VectorMol, mol2::VectorMol)
    # TODO: need bond attribute matching?
    return (b1, b2) -> true
end

function bondmatch(mol::VectorMol, querymol::QueryMol)
    return function (b, qb)
        q = getbond(querymol, qb).query
        return querymatchtree(q, mol, b)
    end
end

bond_match(query::QueryMol, mol::VectorMol) = bond_match(mol, query)


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
        subq = parse(SMARTS, query.second)
        return match_groupquery(mol, subq, i)
    else
        if query.first == :RingSize
            # TODO:
            return query.second in mol.v[query.first][i]
        else
            return mol.v[query.first][i] == query.second
        end
    end
end


function fast_identity_filter(mol1, mol2)
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


function fast_substr_filter(mol1, mol2)
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
