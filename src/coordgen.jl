#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    coordgen!

using Libdl
using Cxx
using coordgenlibs_jll

function coordgen!(mol::GraphMol)
    Libdl.dlopen(libcoordgen, Libdl.RTLD_GLOBAL)
    adir = joinpath(coordgenlibs_jll.artifact_dir, "include/coordgen")
    addHeaderDir(adir)
    # TODO: pass vectors to cpp function
    # @cxx hoge(convert(cxxt"std::vector<int32_t>", [1,2,3,4,5]))
    cxxinclude("sketcherMinimizer.h")
    min_mol = @cxxnew sketcherMinimizerMolecule()
    atoms = []
    bonds = []
    # Atoms
    for nattr in nodeattrs(mol)
        atom = @cxx min_mol->addNewAtom()
        @cxx atom->setAtomicNumber(atomnumber(nattr))
        push!(atoms, atom)
    end
    # Bonds
    for (i, (u, v)) in enumerate(edgesiter(mol))
        bond = @cxx min_mol->addNewBond(atoms[u], atoms[v])
        @cxx bond->setBondOrder(edgeattr(mol, i).order)
        push!(bonds, bond)
    end
    catoms = @cxx min_mol->getAtoms()
    cbonds = @cxx min_mol->getBonds()
    @cxx min_mol->assignBondsAndNeighbors(catoms, cbonds)
    # Stereocenter
    for i in 1:nodecount(mol)
        nattr = nodeattr(mol, i)
        nattr.stereo in (:clockwise, :anticlockwise) || continue
        direction = nattr.stereo === :clockwise
        if degree(mol, i) == 4
            f, s, t, _ = sort(collect(adjacencies(mol, i)))
        elseif degree(mol, i) == 3
            f, s, t = sort(collect(adjacencies(mol, i)))
            if i < f || (i > s && i < t) # implicit H is the 1st or 3rd
                direction = !direction
            end
        else
            continue
        end
        cxx"""
            sketcherMinimizerAtomChiralityInfo getchiralityinfo(
                sketcherMinimizerAtom* lookingFrom,
                sketcherMinimizerAtom* atom1,
                sketcherMinimizerAtom* atom2,
                bool direction
            ) {
                sketcherMinimizerAtomChiralityInfo::sketcherMinimizerChirality d = direction ? sketcherMinimizerAtomChiralityInfo::clockwise : sketcherMinimizerAtomChiralityInfo::counterClockwise;
                return sketcherMinimizerAtomChiralityInfo {
                    lookingFrom, atom1, atom2, d};
            }
        """
        ci = @cxx getchiralityinfo(atoms[f], atoms[s], atoms[t], direction)
        @cxx atoms[i]->setStereoChemistry(ci)
        @cxx atoms[i]->setAbsoluteStereoFromChiralityInfo()
    end
    # StereoBond
    for (i, (u, v)) in enumerate(edgesiter(mol))
        eattr = edgeattr(mol, i)
        eattr.stereo in (:cis, :trans) || continue
        iscis = eattr.stereo === :cis
        fs = Int[]
        for n in (u, v)
            incs = [inc for inc in incidences(mol, n) if inc != i]
            if length(incs) == 1 && incs[1] > i  # lower indexed implicit H
                iscis = !iscis
            end
            push!(fs, sort(incs)[1])
        end
        cxx"""
            sketcherMinimizerBondStereoInfo getstereoinfo(
                sketcherMinimizerAtom* atom1,
                sketcherMinimizerAtom* atom2,
                bool iscis
            ) {
                sketcherMinimizerBondStereoInfo::sketcherMinimizerBondStereo c = iscis ? sketcherMinimizerBondStereoInfo::cis : sketcherMinimizerBondStereoInfo::trans;
                return sketcherMinimizerBondStereoInfo {atom1, atom2, c};
            }
        """
        si = @cxx getstereoinfo(atoms[fs[1]], atoms[fs[2]], iscis)
        @cxx bonds[i]->setStereoChemistry(si)
        @cxx bonds[i]->setAbsoluteStereoFromStereoInfo()
    end
    # Optimize
    minimizer = @cxxnew sketcherMinimizer()
    @cxx minimizer->initialize(min_mol)
    @cxx minimizer->runGenerateCoordinates()
    # Output
    coords = zeros(Float64, nodecount(mol), 2)
    for i in 1:nodecount(mol)
        p = @cxx atoms[i]->getCoordinates()
        px = @cxx p->x()
        py = @cxx p->y()
        coords[i, :] = [px, py]
    end
    notation = Int[]
    for i in 1:edgecount(mol)
        hasstereo = @cxx bonds[i]->hasStereochemistryDisplay
        if !hasstereo
            push!(notation, 0)
            continue
        end
        iswedge = @cxx bonds[i]->isWedge
        isrev = @cxx bonds[i]->isReversed
        push!(notation, iswedge ? 1 : 6)
        u, v = mol.edges[i]
        if isrev
            mol.edges[i] = (v, u)
        end
    end
    mol.cache[:coords2d] = coords
    mol.cache[:coordgenbond] = notation
end
