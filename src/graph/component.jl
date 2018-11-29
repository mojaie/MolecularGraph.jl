#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    connected_components,
    two_edge_connected


function dfswalk(graph, v, visited)
    push!(visited, v)
    for nbr in neighborkeys(graph, v)
        if !(nbr in visited)
            dfswalk(graph, nbr, visited)
        end
    end
end


function connected_components(graph::AbstractUGraph)
    components = Set{Int}[]
    nodeset = nodekeys(graph)
    visited = Set{Int}()
    while !isempty(nodeset)
        v = pop!(nodeset)
        dfswalk(graph, v, visited)
        push!(components, copy(visited))
        setdiff!(nodeset, visited)
        empty!(visited)
    end
    return components
end


function two_edge_connected(graph::AbstractUGraph)
    brs = Int[]
    for conn in connected_components(graph)
        subg = inducedsubgraph(graph, conn)
        root = pop!(nodekeys(subg))
        append!(brs, bridges(subg, root))
    end
    bicomp = GUGraphView(graph, nodekeys(graph), setdiff(edgekeys(graph), brs))
    return Set{Int}[c for c in connected_components(bicomp) if length(c) > 1]
end
