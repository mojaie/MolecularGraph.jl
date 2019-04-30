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

- connected(Bool): if true, apply connected MCS constraint.
- topological(Bool): if true, apply topological constraint.
- diameter(Int): distance cutoff for topological constraint.
- tolerance(Int): distance mismatch tolerance for topological constraint.
- timeout(Int): abort calculation and return suboptimal results if the execution
time has reached the given value (default=60, in seconds).
- targetsize(Int): abort calculation and return suboptimal result so far if the
given mcs size achieved.
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
