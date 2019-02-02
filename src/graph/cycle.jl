#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    cycleedges,
    minimumcycles


function findcycleedge(graph, v, pred, cy)
    for (nbr, edge) in neighbors(graph, v)
        if !(nbr in keys(pred))
            # New node
            pred[nbr] = v
            findcycleedge(graph, nbr, pred, cy)
        elseif nbr != pred[v] && !(edge in cy)
            # Cycle found
            push!(cy, edge)
        end
    end
end


function cycleedges(graph::UDGraph, root)
    pred = Dict{Int, Union{Int,Nothing}}(root => nothing)
    cy = Int[]
    findcycleedge(graph, root, pred, cy)
    return cy
end


function minimumcycles(graph::UDGraph)
    mincycs = Vector{Int}[]
    for biconn in two_edge_connected(graph)
        subg = nodesubgraph(graph, biconn)
        for cyc in mincyclebasis(subg)
            push!(mincycs, cyc)
        end
    end
    return mincycs
end


function mincyclebasis(graph::UDGraph)
    # de Pina algorithm re-interpreted by Kavitha et al.
    cycles = Vector{Int}[]
    root = pop!(nodekeys(graph))
    S = [Set(e) for e in cycleedges(graph, root)]
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
    gnodes = nodekeys(graph)
    G = UDSubgraph(graph, gnodes, setdiff(edgekeys(graph), S))
    U, nmap, emap = shallowmerge(G, G)
    nrev = Dict(v => k for (k, v) in nmap)
    for s in S
        e = getedge(graph, s)
        e1 = connect(e, e.u, nmap[e.v])
        e2 = connect(e, e.v, nmap[e.u])
        maxe = maximum(edgekeys(U))
        updateedge!(U, e1, maxe + 1)
        updateedge!(U, e2, maxe + 2)
    end
    ps = []
    pls = []
    for n in gnodes
        sp = shortestpath(U, n, nmap[n])
        push!(ps, sp)
        push!(pls, length(sp))
    end
    minpath = ps[argmin(pls)]
    ns = Int[n in gnodes ? n : nrev[n] for n in minpath]
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
