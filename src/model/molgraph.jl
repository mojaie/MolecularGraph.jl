#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    AbstractMolecule,
    Annotation,
    Molecule,
    MutableMolecule,
    nullmol,
    getatom,
    getbond,
    neighbors,
    neighborcount,
    atomcount,
    bondcount,
    updateatom!,
    updatebond!,
    required_annotation

import ..GraphModel: neighbors, neighborcount


abstract type AbstractMolecule end
abstract type Annotation end


struct MutableMolecule <: AbstractMolecule
    graph::MutableUDGraph{Atom,Bond}
    attribute::Dict

    function MutableMolecule()
        new(MutableUDGraph{Atom,Bond}(), Dict())
    end
end


struct Molecule <: AbstractMolecule
    graph::UDGraph{Atom,Bond}
    v::Dict{Symbol, Array}
    annotation::Dict{Symbol, Annotation}
    attribute::Dict
end

function Molecule(mol::MutableMolecule)
    Molecule(UDGraph(mol.graph), Dict(), Dict(), mol.attribute)
end

function Molecule(nodes::Vector{Atom}, edges::Vector{Bond})
    # do not use `fill`
    adj = [Dict() for i in 1:length(nodes)]
    for (i, e) in enumerate(edges)
        adj[e.u][e.v] = i
        adj[e.v][e.u] = i
    end
    graph = UDGraph{Atom,Bond}(nodes, edges, adj)
    Molecule(graph, Dict(), Dict(), Dict())
end


nullmol() = Molecule(Vector{Atom}[], Vector{Bond}[])


getatom(mol::AbstractMolecule, idx) = getnode(mol.graph, idx)

getbond(mol::AbstractMolecule, u, v) = getedge(mol.graph, u, v)

neighbors(mol::AbstractMolecule, idx) = neighbors(mol.graph, idx)

atomcount(mol::AbstractMolecule) = nodecount(mol.graph)
bondcount(mol::AbstractMolecule) = edgecount(mol.graph)

neighborcount(mol::AbstractMolecule, idx) = length(neighbors(mol.graph, idx))


function updateatom!(mol::MutableMolecule, atom, idx)
    updatenode!(mol.graph, atom, idx)
end


function updatebond!(mol::MutableMolecule, bond, idx)
    updateedge!(mol.graph, bond, idx)
end

updatebond!(m::MutableMolecule, bond, u, v) = updatedge!(m.graph, bond, u, v)


function required_annotation(mol::Molecule, annot)
    if !(annot in keys(mol.annotation))
        throw(ErrorException("$(annot) is not available"))
    end
end
