#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


"""
    lgvmatch(g, h, grev, hrev, vmatch, ematch) -> Function

A line graph vmatch function for MCES calculation.
"""
function lgvmatch(
        g::SimpleGraph{T}, h::SimpleGraph{T}, grev::Dict{T,Edge{T}}, hrev::Dict{T,Edge{T}},
        vmatch::F1, ematch::F2) where {T,F1,F2}
    return function (lgn, lhn)
        ge = grev[lgn]
        he = hrev[lhn]
        ematch(ge, he) || return false
        return (
            (vmatch(src(ge), src(he)) && vmatch(dst(ge), dst(he)))
            || (vmatch(src(ge), dst(he)) && vmatch(dst(ge), src(he)))
        )
    end
end


"""
    lgematch(g, h, gsh, hsh, vmatch) -> Function

A line graph ematch function for MCES calculation.
"""
function lgematch(
        g::SimpleGraph{T}, h::SimpleGraph{T},
        gsh::Dict{Edge{T},T}, hsh::Dict{Edge{T},T}, vmatch::F) where {T,F}
    return (lge, lhe) -> vmatch(gsh[lge], hsh[lhe])
end



function lgvmatchvecgen(
        revmap::Dict{T,Edge{T}}, vmatchvec::F1, ematchvec::F2) where {T,F1,F2}
    return function (i)
        e = revmap[i]
        u, v = (src(e), dst(e))
        uvec = vmatchvec(u)
        vvec = vmatchvec(v)
        evec = ematchvec(u, v)
        # uvec 7 bit + vvec 7 bit + evec 2 bit
        return (uvec <= vvec ? evec * 14400 + (uvec * 120 + vvec)
            : evec * 14400 + (vvec * 120 + uvec))
    end
end


function lgematchvecgen(
        lg::SimpleGraph{T}, shared::Dict{Edge{T},T}, vmatchvec::F) where {T,F}
    return function (u, v)
        has_edge(lg, u, v) || return 0
        return vmatchvec(shared[u_edge(lg, u, v)])
    end
end




"""
    delta_edges(g::SimpleGraph{T}, vmatch=(n1,n2)->true) -> Vector{Tuple{Edge{T},Edge{T},Edge{T}}}

Returns 3-tuples of edges that forms a triangle with uniform vertex attributes.

Used for Δ-Y transformation check in edge-induced subgraph isomorphism.
"""
function delta_edges(g::SimpleGraph{T}, vmatch=(n1,n2)->true) where T
    ts = Set{Tuple{Edge{T},Edge{T},Edge{T}}}()
    for n in vertices(g)
        nbrs = neighbors(g, n)
        length(nbrs) < 2 && continue
        for (i, j) in combinations(length(nbrs)) # find triangles
            u, v = (nbrs[i], nbrs[j])
            if has_edge(g, u, v) && vmatch(n, u) && vmatch(n, v) && vmatch(u, v)
                push!(ts, tuple(u_edge(g, n, u), u_edge(g, n, v), u_edge(g, u, v)))
            end
        end
    end
    return collect(ts)
end


"""
    y_edges(g::SimpleGraph{T}, vmatch=(n1,n2)->true) -> Vector{Tuple{Edge{T},Edge{T},Edge{T}}}

Returns 3-tuples of edges that forms a star graph with uniform vertex attributes.

Used for Δ-Y transformation check in edge-induced subgraph isomorphism.
"""
function y_edges(g::SimpleGraph{T}, vmatch=(n1,n2)->true) where T
    ts = Set{Tuple{Edge{T},Edge{T},Edge{T}}}()
    for n in vertices(g)
        nbrs = neighbors(g, n)
        length(nbrs) < 3 && continue
        for (i, j, k) in combinations(length(nbrs), 3)
            u, v, w = (nbrs[i], nbrs[j], nbrs[k])
            if vmatch(n, u) && vmatch(n, v) && vmatch(n, w)
                push!(ts, tuple(u_edge(g, n, u), u_edge(g, n, v), u_edge(g, n, w)))
            end
        end
    end
    return collect(ts)
end


"""
    delta_y_test(mapping, gtri, gy, htri, hy) -> Bool

Returns whether the edge-induced isomorphism mapping does not have Δ-Y mismatches.
"""
function delta_y_test(
        mapping::Dict{Edge{T},Edge{T}},  # edge mapping g(e) => h(e)
        gtri::Vector{Tuple{Edge{T},Edge{T},Edge{T}}},
        gy::Vector{Tuple{Edge{T},Edge{T},Edge{T}}},
        htri::Vector{Tuple{Edge{T},Edge{T},Edge{T}}},
        hy::Vector{Tuple{Edge{T},Edge{T},Edge{T}}}) where T
    revmap = Dict(v => k for (k, v) in mapping)
    for t in gtri  # triangle edges
        issubset(t, keys(mapping)) || continue
        for y in hy
            issubset(y, keys(revmap)) || continue
            if issetequal(t, [revmap[e] for e in y])
                return false
            end
        end
    end
    for t in htri
        issubset(t, keys(revmap)) || continue
        for y in gy
            issubset(y, keys(mapping)) || continue
            if issetequal(t, [mapping[e] for e in y])
                return false
            end
        end
    end
    return true
end


"""
    delta_y_correction!(mapping, gtri, gy, htri, hy) -> Nothing

Remove edges from the edge-induced isomorphism mapping to correct Δ-Y mismatches.

If vertex attribute matching is enabled, any attribute mismatch would break isomorphism of triangles,
so this function only deals with triangles and stars that have uniform vertex attributes.
Note that any edge attribute mismatch would not affect isomorphism of other two edges in the triplets.
"""
function delta_y_correction!(
        mapping::Dict{Edge{T},Edge{T}},  # edge mapping g(e) => h(e)
        gtri::Vector{Tuple{Edge{T},Edge{T},Edge{T}}},
        gy::Vector{Tuple{Edge{T},Edge{T},Edge{T}}},
        htri::Vector{Tuple{Edge{T},Edge{T},Edge{T}}},
        hy::Vector{Tuple{Edge{T},Edge{T},Edge{T}}}) where T
    revmap = Dict(v => k for (k, v) in mapping)
    for t in gtri  # triangle edges
        issubset(t, keys(mapping)) || continue
        for y in hy
            issubset(y, keys(revmap)) || continue
            if issetequal(t, [revmap[e] for e in y])
                delete!(mapping, t[1]) # remove arbitrary one edge in the triangle
            end
        end
    end
    for t in htri
        issubset(t, keys(revmap)) || continue
        for y in gy
            issubset(y, keys(mapping)) || continue
            if issetequal(t, [mapping[e] for e in y])
                delete!(mapping, revmap[t[1]])
            end
        end
    end
end
