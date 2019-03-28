#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    mincycles, circuitrank,
    node_cyclemem, edge_cyclemem


function find_cotree_edge(graph, v, pred, cy)
    for (nbr, edge) in neighbors(graph, v)
        if !(nbr in keys(pred))
            # New node
            pred[nbr] = v
            find_cotree_edge(graph, nbr, pred, cy)
        elseif nbr != pred[v] && !(edge in cy)
            # Cycle found
            push!(cy, edge)
        end
    end
end


function cotree_edges(graph::UndirectedGraph, root)
    pred = Dict{Int, Union{Int,Nothing}}(root => nothing)
    cy = Int[]
    find_cotree_edge(graph, root, pred, cy)
    return cy
end


"""
    mincycles(graph::UndirectedGraph) -> Vector{Vector{Int}}

Calculate minimum cycle basis (also known as Smallest Set of Smallest Rings
in the context of molecular graph theory).
"""
@cache function mincycles(graph::UndirectedGraph)
    mincycs = Vector{Int}[]
    for biconn in two_edge_connected(graph)
        subg = nodesubgraph(graph, biconn)
        for cyc in mincyclebasis(subg)
            push!(mincycs, cyc)
        end
    end
    return mincycs
end


circuitrank(graph::UndirectedGraph) = length(mincycles(graph))


@cache function node_cyclemem(graph::UndirectedGraph)
    nmap = Dict(n => Int[] for n in nodekeys(graph))
    for (i, cyc) in enumerate(mincycles(graph))
        for n in cyc
            push!(nmap[n], i)
        end
    end
    nodes = [nmap[n] for n in nodekeys(graph)]
    return nodes
end


@cache function edge_cyclemem(graph::UndirectedGraph)
    emap = Dict(e => Int[] for e in edgekeys(graph))
    for (i, cyc) in enumerate(mincycles(graph))
        for e in edgekeys(nodesubgraph(graph, cyc))
            push!(emap[e], i)
        end
    end
    edges = [emap[e] for e in edgekeys(graph)]
    return edges
end


function mincyclebasis(graph::UndirectedGraph)
    # de Pina algorithm re-interpreted by Kavitha et al.
    cycles = Vector{Int}[]
    root = pop!(nodeset(graph))
    S = [Set(e) for e in cotree_edges(graph, root)]
    N = length(S)
    for k in 1:N
        minnodes, minedges = findmincycle(graph, S[k])
        push!(cycles, minnodes)
        for i in (k + 1):N
            if length(intersect(S[i], minedges)) % 2 == 1
                S[i] = symdiff(S[i], S[k])
            end
        end
    end
    return cycles
end


function findmincycle(graph, S)
    G = mapgraph(edgesubgraph(graph, setdiff(edgeset(graph), S)))
    U, nmap, emap = shallowmerge(G, G)
    nrev = Dict(v => k for (k, v) in nmap)
    for s in S
        e = getedge(graph, s)
        e1 = setnodes(e, e.u, nmap[e.v])
        e2 = setnodes(e, nmap[e.u], e.v)
        maxe = maximum(edgekeys(U))
        updateedge!(U, e1, maxe + 1)
        updateedge!(U, e2, maxe + 2)
    end
    ps = []
    pls = []
    for n in nodeset(graph)
        sp = shortestpath(U, n, nmap[n])
        push!(ps, sp)
        push!(pls, length(sp))
    end
    minpath = ps[argmin(pls)]
    ns = Int[n in nodeset(graph) ? n : nrev[n] for n in minpath]
    es = [neighbors(graph, ns[n])[ns[n+1]] for n in 1:(length(ns) - 1)]
    pop!(ns)
    return canonicalize(ns), es
end


function canonicalize(nodes)
    """Align cycle indices to start from lowest index and following one of
    neighbors that have the lower index
    """
    (fst, fstidx) = findmin(nodes)
    succidx = fstidx == lastindex(nodes) ? 1 : fstidx + 1
    succ = nodes[succidx]
    predidx = fstidx == 1 ? lastindex(nodes) : fstidx - 1
    pred = nodes[predidx]
    cp = succ < pred ? copy(nodes) : reverse(copy(nodes))
    while cp[1] != fst
        pushfirst!(cp, pop!(cp))
    end
    return cp
end
