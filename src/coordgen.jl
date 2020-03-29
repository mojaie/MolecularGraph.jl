#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    coordgen

using Libdl
using Cxx
using coordgenlibs_jll

function coordgen(mol::GraphMol)
    Libdl.dlopen(libcoordgen, Libdl.RTLD_GLOBAL)
    adir = joinpath(coordgenlibs_jll.artifact_dir, "include/coordgen")
    addHeaderDir(adir)
    cxxinclude("sketcherMinimizer.h")
    min_mol = @cxxnew sketcherMinimizerMolecule()
    atoms = []
    for node in nodeattrs(mol)
        a = @cxx min_mol->addNewAtom()
        @cxx a->setAtomicNumber(atomnumber(node))
        push!(atoms, a)
    end
    order = bondorder(mol)
    for (i, (u, v)) in enumerate(edgesiter(mol))
        b = @cxx min_mol->addNewBond(atoms[u], atoms[v])
        @cxx b->setBondOrder(order[i])
    end
    minimizer = @cxxnew sketcherMinimizer()
    @cxx minimizer->initialize(min_mol)
    @cxx minimizer->runGenerateCoordinates()
    coords = zeros(Float64, nodecount(mol), 2)
    for i in 1:nodecount(mol)
        p = @cxx atoms[i]->getCoordinates()
        px = @cxx p->x()
        py = @cxx p->y()
        coords[i, :] = [px, py]
    end
    return coords
end
