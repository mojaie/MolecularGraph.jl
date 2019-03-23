#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    atomsubstr,
    bondsubstr

# TODO: carry over original annotations

function atomsubstr(mol::MapMol, atoms)
    subg = vectorgraph(nodesubgraph(mol, atoms))
    return mapmol(subg.nodes, subg.edges)
end

function atomsubstr(mol::VectorMol, atoms)
    subg = nodesubgraph(mol, atoms)
    onodes = sort(nodekeys(subg))
    oedges = sort(edgekeys(subg))
    subg = vectorgraph(subg)
    newmol = vectormol(subg.nodes, subg.edges)
    for (k, vec) in mol.vector
        # TODO: workaround
        keys = occursin("Bond", string(k)) ? oedges : onodes
        newmol.vector[k] = vec[keys]
    end
    return newmol
end


function bondsubstr(mol::MapMol, bonds)
    subg = vectorgraph(edgesubgraph(mol, bonds))
    return mapmol(subg.nodes, subg.edges)
end

function bondsubstr(mol::VectorMol, bonds)
    subg = edgesubgraph(mol, bonds)
    onodes = sort(nodekeys(subg))
    oedges = sort(edgekeys(subg))
    subg = vectorgraph(subg)
    newmol = vectormol(subg.nodes, subg.edges)
    for (k, vec) in mol.vector
        # TODO: workaround
        keys = occursin("Bond", string(k)) ? oedges : onodes
        newmol.vector[k] = vec[keys]
    end
    return newmol
end
