#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    shallowmerge


function shallowmerge(G::UDGraph, H::UDGraph)
    U = similarmap(G)
    for (i, n) in nodesiter(G)
        U.nodes[i] = n
        U.adjacency[i] = Dict()
    end
    for (i, e) in edgesiter(G)
        U.edges[i] = e
        U.adjacency[e.u][e.v] = i
        U.adjacency[e.v][e.u] = i
    end
    nmax = maximum(nodekeys(U))
    nmap = Dict(i => i + nmax for i in nodekeys(H))
    emax = maximum(edgekeys(U))
    emap = Dict(i => i + emax for i in edgekeys(H))
    for (i, n) in nodesiter(H)
        U.nodes[nmap[i]] = n
        U.adjacency[nmap[i]] = Dict()
    end
    for (i, e) in edgesiter(H)
        U.edges[emap[i]] = e
        U.adjacency[nmap[e.u]][nmap[e.v]] = emap[i]
        U.adjacency[nmap[e.v]][nmap[e.u]] = emap[i]
    end
    return U, nmap, emap
end
