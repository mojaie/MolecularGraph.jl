#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    VectorMol, QueryMol, MapMol,
    vectormol, querymol, mapmol,
    SubstructureView,
    SDFile, SMILES, SMARTS,
    getatom, getbond, atomcount, bondcount,
    updateatom!, updatebond!, unlinkatoms, unlinkbonds


struct VectorMol{A<:Atom,B<:Bond} <: GeneralMol
    graph::VectorGraph{A,B}
    attribute::Dict{Symbol,Any}
    cache::Dict{Symbol,Any}

    function VectorMol{A,B}(graph::VectorGraph{A,B}) where {A<:Atom,B<:Bond}
        new(graph, Dict(), Dict())
    end
end

VectorMol{A,B}() where {A<:Atom,B<:Bond} = VectorMol{A,B}(vectorgraph(A,B))


"""
    vectormol(::Type{A}, ::Type{B}) where {A<:Atom,B<:Bond} -> VectorMol{A,B}

Generate empty `VectorMol` that has atoms and bonds with the given types.
"""
vectormol(::Type{A}, ::Type{B}) where {A<:Atom,B<:Bond} = VectorMol{A,B}()

"""
    vectormol(atoms::Vector{A}, bonds::Vector{B}) -> VectorMol{A,B}

Generate `VectorMol` that has the given atom objects and edge objects.
"""
function vectormol(nodes::Vector{A}, edges::Vector{B}) where {A<:Atom,B<:Bond}
    mol = VectorMol{A,B}()
    for (i, node) in enumerate(nodes)
        push!(mol.graph.nodes, node)
        push!(mol.graph.neighbormap, Dict())
    end
    for (i, edge) in enumerate(edges)
        push!(mol.graph.edges, edge)
        mol.graph.neighbormap[edge.u][edge.v] = i
        mol.graph.neighbormap[edge.v][edge.u] = i
    end
    return mol
end

"""
    vectormol(mol::GraphMol; clone=false) -> VectorMol

Convert the given molecule into a new `VectorMol`. See [`vectorgraph`](@ref)
for the details.
"""
function vectormol(mol::GeneralMol)
    A = nodetype(mol)
    B = edgetype(mol)
    newmol = VectorMol{A,B}(vectorgraph(mol))
    merge!(newmol.attribute, mol.attribute)
    return newmol
end



struct QueryMol{A<:QueryAtom,B<:QueryBond} <: GraphMol
    graph::MapGraph{A,B}
    connectivity::Array{Array{Int}}

    function QueryMol{A,B}() where {A<:QueryAtom,B<:QueryBond}
        new(mapgraph(A,B), [])
    end
end

"""
    querymol(::Type{A}, ::Type{B}
        ) where {A<:QueryAtom,B<:QueryBond} -> QueryMol{N,E}

Generate empty `QueryMol` that has atoms and bonds with the given types.
"""
querymol(::Type{A}, ::Type{B}
    ) where {A<:QueryAtom,B<:QueryBond} = QueryMol{A,B}()



struct MapMol{A<:Atom,B<:Bond} <: GraphMol
    # TODO: deprecated
    graph::MapGraph{A,B}
    attribute::Dict{Symbol,String}

    function MapMol{A,B}(graph::MapGraph{A,B}) where {A<:Atom,B<:Bond}
        new(graph, Dict())
    end
end

MapMol{A,B}() where {A<:Atom,B<:Bond} = MapMol{A,B}(mapgraph(A,B))

mapmol(::Type{A}, ::Type{B}) where {A<:Atom,B<:Bond} = MapMol{A,B}()

function mapmol(atoms::Vector{A}, bonds::Vector{B}) where {A<:Atom,B<:Bond}
    mol = MapMol{A,B}()
    for (i, atom) in enumerate(atoms)
        mol.graph.nodes[i] = atom
        mol.graph.neighbormap[i] = Dict()
    end
    for (i, bond) in enumerate(bonds)
        mol.graph.edges[i] = bond
        mol.graph.neighbormap[bond.u][bond.v] = i
        mol.graph.neighbormap[bond.v][bond.u] = i
    end
    return mol
end

function mapmol(mol::GraphMol)
    A = nodetype(mol)
    B = edgetype(mol)
    newmol = MapMol{A,B}(mapgraph(mol))
    merge!(newmol.attribute, mol.attribute)
    return newmol
end



struct SubstructureView{T<:UndirectedGraph} <: GeneralMolView
    graph::SubgraphView{T}
    attribute::Dict{Symbol,Any}
    cache::Dict{Symbol,Any}
end


function MolecularGraphModel.nodesubgraph(mol::GeneralMol, nodes)
    subg = nodesubgraph(mol.graph, nodes)
    return SubstructureView(subg, mol.attribute, mol.cache)
end


function MolecularGraphModel.edgesubgraph(mol::GeneralMol, edges)
    subg = edgesubgraph(mol.graph, edges)
    return SubstructureView(subg, mol.attribute, mol.cache)
end


Base.getindex(mol::GraphMol, sym::Symbol) = eval(Expr(:call, sym, mol))
Base.getindex(
    mol::GraphMol, k1::Symbol, k2::Symbol, K::Symbol...
) = hcat(eval(Expr(:call, sym, mol)) for k in [k1, k2, K...])


# Aliases

SDFile = VectorMol{SDFileAtom,SDFileBond}
SMILES = VectorMol{SmilesAtom,SmilesBond}
SMARTS = QueryMol{SmartsAtom,SmartsBond}

getatom = getnode
getbond = getedge
atomcount = nodecount
bondcount = edgecount
updateatom! = updatenode!
updatebond! = updateedge!
unlinkatoms = unlinknodes
unlinkbonds = unlinkedges
