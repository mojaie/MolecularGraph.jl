#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    shallowmerge


function shallowmerge(G::UndirectedGraph, H::UndirectedGraph)
    U = mapgraph(G)
    for (i, n) in nodesiter(G)
        U.nodes[i] = n
        U.neighbormap[i] = Dict()
    end
    for (i, e) in edgesiter(G)
        U.edges[i] = e
        U.neighbormap[e.u][e.v] = i
        U.neighbormap[e.v][e.u] = i
    end
    nmax = maximum(nodekeys(U))
    nmap = Dict(i => i + nmax for i in nodekeys(H))
    emax = maximum(edgekeys(U))
    emap = Dict(i => i + emax for i in edgekeys(H))
    for (i, n) in nodesiter(H)
        U.nodes[nmap[i]] = n
        U.neighbormap[nmap[i]] = Dict()
    end
    for (i, e) in edgesiter(H)
        newe = setnodes(e, nmap[e.u], nmap[e.v])
        eidx = emap[i]
        U.edges[eidx] = newe
        U.neighbormap[newe.u][newe.v] = eidx
        U.neighbormap[newe.v][newe.u] = eidx
    end
    return U, nmap, emap
end
