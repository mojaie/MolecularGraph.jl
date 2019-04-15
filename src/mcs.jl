#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    mcismol, mcismolsize,
    mcesmol, mcesmolsize


"""
    mcismol(mol1::UndirectedGraph, mol2::UndirectedGraph; kwargs...
        ) -> MCSResults

Compute maximum common induced substructure (MCIS) of mol1 and mol2.

## Keyword arguments

- timeout(Int): abort calculation and return suboptimal results if the execution
time exceeded the value (default=60, in seconds)
- c_clique_constraint(Bool): if true, calculate connected MCS.
"""
function mcismol(mol1::UndirectedGraph, mol2::UndirectedGraph; kwargs...)
    afunc = atommatch(mol1, mol2)
    bfunc = bondmatch(mol1, mol2)
    return findmcis(
        mol1, mol2, nodematcher=afunc, edgematcher=bfunc; kwargs...)
end

mcismolsize(mol1, mol2; kwargs...
    ) = length(mcismol(mol1, mol2; kwargs...).mapping)


"""
    mcesmol(mol1::UndirectedGraph, mol2::UndirectedGraph; kwargs...
        ) -> MCSResults

Compute maximum common edge induced substructure (MCES) of mol1 and mol2.
"""
function mcesmol(mol1::UndirectedGraph, mol2::UndirectedGraph; kwargs...)
    afunc = atommatch(mol1, mol2)
    bfunc = bondmatch(mol1, mol2)
    return findmces(
        mol1, mol2, nodematcher=afunc, edgematcher=bfunc; kwargs...)
end

mcesmolsize(mol1, mol2; kwargs...
    ) = length(mcesmol(mol1, mol2; kwargs...).mapping)

# TODO: subgraphview
