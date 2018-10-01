module UndirectedGraph

export Graph, node, edge, updatenode!, updateedge!


struct Graph
    adjacency::Dict
    nodes::List
    function Graph(nodes::List, edges::List)
        new(adjacency, nodes)
    end
end


function node(graph::Graph, idx)
    graph.nodes[idx]
end


function edge(graph::Graph, u, v)
    graph.adjacency[u][v]
end


function updatenode!(graph::Graph, idx::Int32, attr::Dict)
    graph.nodes[idx] = attr
end


function updateedge!(graph::Graph, u::Int32, v::Int32, attr::Dict)
    graph.adjacency[u][v] = attr
    graph.adjacency[v][u] = attr
end


end  # module
