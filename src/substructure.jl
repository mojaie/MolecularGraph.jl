#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    STRUCT_MATCH_SETTING,
    is_identical,
    is_substruct,
    is_superstruct,
    substruct_mappings,
    struct_match,
    preprocess,
    fast_identity_filter,
    fast_substr_filter,
    deltaYmap!


const STRUCT_MATCH_SETTING = Dict(
    :prefilter => (a, b) -> true,
    :preprocess => identity,
    :postprocess => identity,
    :vf2conf => Dict(
        :mode => :isomorphic,
        :depthlimit => 100,
        :node_match => nothing,
        :edge_match => nothing
    )
)


function is_identical(mol1, mol2)
    setting = copy(STRUCT_MATCH_SETTING)
    setting[:prefilter] = fast_identity_filter
    setting[:preprocess] = preprocess
    iterate(struct_match(mol1, mol2, setting)) !== nothing
end


function is_substruct(mol1, mol2)
    setting = copy(STRUCT_MATCH_SETTING)
    setting[:vf2conf][:mode] = :subgraph
    setting[:prefilter] = fast_substr_filter
    setting[:preprocess] = preprocess
    iterate(struct_match(mol2, mol1, setting)) !== nothing
end

is_superstruct(mol1, mol2) = is_substruct(mol2, mol1)


function substruct_mappings(mol1, mol2)
    setting = copy(STRUCT_MATCH_SETTING)
    setting[:vf2conf][:mode] = :subgraph
    setting[:prefilter] = fast_substr_filter
    setting[:preprocess] = preprocess
    struct_match(mol2, mol2, setting)
end


function struct_match(mol1, mol2, setting)
    Channel() do channel
        if !setting[:prefilter](mol1, mol2)
            return
        end
        lg1 = setting[:preprocess](mol1)
        lg2 = setting[:preprocess](mol2)
        vf2conf = copy(setting[:vf2conf])
        vf2conf[:node_match] = node_match(mol1, mol2, lg1, lg2)
        vf2conf[:edge_match] = edge_match(mol1, mol2, lg1, lg2)
        state = VF2State(lg1, lg2, vf2conf)
        for mapping in Channel(c -> vf2match!(c, state, nothing, nothing))
            # Delta-Y transformation check
            dy1 = deltaYmap(mapping, mol1, mol2)
            dy2 = deltaYmap(Dict(v => k for (k, v) in mapping), mol2, mol1)
            for m in keys(dy1)
                delete!(mapping, m[1])
            end
            for m in values(dy2)
                delete!(mapping, m[1])
            end
            put!(channel, mapping)
        end
    end
end


function preprocess(mol)
    # remove and annotate Hs and trivials
    required_annotation(mol, :Topology)
    required_annotation(mol, :Default)
    linegraph(mol.graph)
end


function node_match(mol1, mol2, lg1, lg2)
    function (g, h)
        sym1 = mol1.v[:Symbol]
        pi1 = mol1.v[:Pi]
        sym2 = mol2.v[:Symbol]
        pi2 = mol2.v[:Pi]
        p = getnode(lg1, g)
        q = getnode(lg2, h)

        ((sym1[p.n1] == sym2[q.n1] && sym1[p.n2] == sym2[q.n2]
         && pi1[p.n1] == pi2[q.n1] && pi1[p.n2] == pi2[q.n2])
        || (sym1[p.n1] == sym2[q.n2] && sym1[p.n2] == sym2[q.n1]
         && pi1[p.n1] == pi2[q.n2] && pi1[p.n2] == pi2[q.n1]))
    end
end


function edge_match(mol1, mol2, lg1, lg2)
    function (g, h)
        sym1 = mol1.v[:Symbol]
        pi1 = mol1.v[:Pi]
        sym2 = mol2.v[:Symbol]
        pi2 = mol2.v[:Pi]
        p = getedge(lg1, g)
        q = getedge(lg2, h)
        sym1[p.common] == sym2[q.common] && pi1[p.common] == pi2[q.common]
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


function deltaYmap(mapping, mol1, mol2)
    # Return edges should be removed
    res = []
    for r in mol1.annotation[:Topology].rings
        if length(r) != 3
            continue
        end
        symr = mol1.v[:Symbol][r]
        pir = mol1.v[:Pi][r]
        if all(symr .== symr[1]) && all(pir .== pir[1])
            cyc = [(1, 2), (2, 3), (3, 1)]
            edges = [neighbors(mol1, r[u])[r[v]] for (u, v) in cyc]
            emap = Dict(e => mapping[e] for e in edges if e in keys(mapping))
            if length(emap) < 3
                continue
            end
            m2n = [(e = getbond(mol2, n); (e.u, e.v)) for n in values(emap)]
            if length(intersect(m2n)) == 4  # star-shaped
                push!(res, emap)
            end
        end
    end
    res
end
