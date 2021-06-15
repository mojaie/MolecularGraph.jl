#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export maximummatching, isperfectmatching


function backtrack(pred, n)
    path = [n]
    x = n
    while pred[x] != x
        x = pred[x]
        pushfirst!(path, x)
    end
    return path
end


function findaugpath(G::UndirectedGraph, matching)
    head = nodeset(G)  # Initial exposed node set
    for m in matching
        (u, v) = getedge(G, m)
        setdiff!(head, u, v)
    end
    @debug "head: $(head)"
    pred = Dict{Int,Int}()  # F tree
    for v in head
        pred[v] = v
    end
    while !isempty(head)
        v = pop!(head)
        @debug "  v: $(v)"
        unmarked = setdiff(edgeset(G), copy(matching))
        for e in intersect(incidences(G, v), unmarked)
            w = neighbors(G, v)[e]
            @debug "    w: $(w)"
            if !(w in keys(pred))
                e2s = intersect(incidences(G, w), matching)  # TODO use only() in Julia 1.4
                @assert length(e2s) == 1
                e2 = pop!(e2s)
                x = neighbors(G, w)[e2]
                @debug "      x: $(x)"
                pred[w] = v
                pred[x] = w
                union!(head, x)
            else
                p1 = backtrack(pred, v)
                p2 = backtrack(pred, w)
                length(p2) % 2 == 1 || continue
                if p1[1] != p2[1]
                    P = vcat(p1, reverse(p2))
                    @debug "      P: $(P) (extend tree)"
                    return P
                end
                # Contraction
                root = intersect(p1, p2)[end]
                bpath = union([root], setdiff(p1, p2), reverse(setdiff(p2, p1)))
                Gcont = edgecontraction(G, bpath)
                b = Gcont.vnode
                @debug "      Recursion: $(Gcont.collapsed) -> $(b)"
                Mcont = setdiff(matching, edgeset(nodesubgraph(G, bpath)))
                P = findaugpath(Gcont, Mcont)
                @assert length(P) % 2 == 0
                pos = findfirst(x -> x == b, P)
                if pos === nothing
                    @debug "      P: $(P) (no blossom in augpath)"
                    return P
                end
                # Lifting
                if pos == 1 || pos == length(P)
                    nbrpos = pos == 1 ? 2 : length(P) - 1
                    inc = findedgekey(Gcont, P[pos], P[nbrpos])
                    dst = neighbors(G, P[nbrpos])[inc]
                    if root == dst
                        P[pos] = root
                        @debug "      P: $(P) (blossom root is exposed)"
                        return P
                    end
                else
                    inedge = findedgekey(Gcont, P[pos], P[pos-1])
                    outedge = findedgekey(Gcont, P[pos], P[pos+1])
                    inlet = neighbors(G, P[pos-1])[inedge]
                    outlet = neighbors(G, P[pos+1])[outedge]
                    @assert inlet == root || outlet == root
                    dst = inlet == root ? outlet : inlet
                end
                # Build blossom graph from edges
                bgedges = Int[]
                push!(bgedges, findedgekey(G, bpath[1], bpath[end]))
                for i in 1:(length(bpath) - 1)
                    push!(bgedges, findedgekey(G, bpath[i], bpath[i + 1]))
                end
                bgraph = edgesubgraph(G, bgedges)
                # Find path in the blossom
                r = shortestpathnodes(bgraph, root, dst)
                if length(r) % 2 == 0
                    # find even length path
                    redges = shortestpathedges(bgraph, root, dst)
                    cobg = edgesubgraph(bgraph, setdiff(bgedges, redges))
                    r = shortestpathnodes(cobg, root, dst)
                end
                baug = pos % 2 == 0 ? reverse(r) : r
                P = union(P[1:pos-1], baug, P[pos+1:end])
                @debug "      P: $(P) (lifted)"
                return P
            end
        end
    end
    return Int[]
end


function findmaxmatch(G::UndirectedGraph, matching)
    P = findaugpath(G, matching)
    @assert length(P) % 2 == 0
    if !isempty(P)
        Pedges = Set{Int}()
        for i in 1:length(P)-1
            push!(Pedges, findedgekey(G, P[i], P[i+1]))
        end
        augment = symdiff(Pedges, matching)
        @assert length(matching) + 1 == length(augment)
        return findmaxmatch(G, augment)
    else
        return matching
    end
end


"""
    maximummatching(G::UndirectedGraph; method=:Blossom) -> Set{Int}

Compute maximum cardinality matching by Edmonds' blossom algorithm and return the set of matched edges.
"""
maximummatching(G::UndirectedGraph; method=:Blossom) = findmaxmatch(G, Set{Int}())


"""
    isperfectmatching(G::UndirectedGraph) -> Bool
    isperfectmatching(G::UndirectedGraph, matching::Set{Int}) -> Bool

Return if the given graph has a perfect matching.
"""
isperfectmatching(G::UndirectedGraph) = length(maximummatching(G)) * 2 == nodecount(G)
isperfectmatching(G::UndirectedGraph, matching::Set{Int}) = length(matching) * 2 == nodecount(G)
