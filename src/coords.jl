#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    coords2d, coords3d, coordgen, coordgen!

using coordgenlibs_jll


function coords2d(mol::SimpleMolGraph)
    has_state(mol, :v_coords2d) && return get_state(mol, :v_coords2d)
    coords = zeros(Float64, nv(mol), 2)
    # TODO: if coords not available, throw error and suggest to use coordgen
    for i in vertices(mol)
        coords[i, :] = get_prop(mol, i, :coords)[1:2]
    end
    return coords
end


function coords3d(mol::SimpleMolGraph)
    has_state(mol, :v_coords3d) && return get_state(mol, :v_coords3d)
    coords = zeros(Float64, nv(mol), 3)
    # TODO: if coords not available, throw error and suggest to use coordgen
    for i in vertices(mol)
        coords[i, :] = get_prop(mol, i, :coords)[1:3]
    end
    return coords
end


"""
    coordgen(mol::MolGraph) -> Tuple{Array{Int,1},Array{Int,1}}

Generate 2D coordinates by using Schrodinger's coordgenlibs.

This will returns a tuple of `coords` and `styles` arrays. `coords` is a size(n, 2) matrix
where n is atom count, which stores 2D coordinates (x, y) of each atoms.
`styles` is a size e vector of wedge notation of stereobond, where e is bond count.
"""
function coordgen(mol::SimpleMolGraph)
    # properties and descriptors
    bondorder_ = bond_order(mol)
    stereocenters = get_prop(mol, :stereocenter)
    stereobonds = get_prop(mol, :stereobond)

    minmol = ccall((:getSketcherMinimizer, libcoordgen), Ptr{Cvoid}, ())
    atoms = Ptr{Cvoid}[]
    bonds = Ptr{Cvoid}[]

    # Atoms
    for a in atom_symbol(mol)
        atom = ccall(
            (:setAtom, libcoordgen), Ptr{Cvoid},
            (Ptr{Cvoid}, Int), minmol, atomnumber(a)
        )
        push!(atoms, atom)
    end

    # Bonds
    for (i, e) in enumerate(edges(mol))
        bond = ccall(
            (:setBond, libcoordgen), Ptr{Cvoid},
            (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Int),
            minmol, atoms[src(e)], atoms[dst(e)], bondorder_[i]
        )
        push!(bonds, bond)
    end

    ccall(
        (:assignBondsAndNeighbors, libcoordgen), Cvoid,
        (Ptr{Cvoid},), minmol
    )

    # Stereocenter
    for (n, stereo) in stereocenters
        ccall(
            (:setStereoCenter, libcoordgen), Cvoid,
            (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Int),
            atoms[n], atoms[stereo[1]], atoms[stereo[2]], atoms[stereo[3]], stereo[4]
        )
    end

    # Stereobond
    for (e, stereo) in stereobonds
        ccall(
            (:setStereoBond, libcoordgen), Cvoid,
            (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Int),
            bonds[edge_rank(mol, e)], atoms[stereo[1]], atoms[stereo[2]], stereo[3]
        )
    end

    # Optimize
    ccall(
        (:runGenerateCoordinates, libcoordgen), Cvoid,
        (Ptr{Cvoid},), minmol
    )

    # Output
    coords = zeros(Float64, nv(mol), 2)
    for i in vertices(mol)
        px = ccall(
            (:getAtomX, libcoordgen), Float32, (Ptr{Cvoid},), atoms[i])
        py = ccall(
            (:getAtomY, libcoordgen), Float32, (Ptr{Cvoid},), atoms[i])
        coords[i, :] = [px, py]
    end
    styles = init_edge_descriptor(Symbol, mol)
    for (i, e) in enumerate(edges(mol))
        hasstereo = ccall(
            (:hasStereochemistryDisplay, libcoordgen), Bool,
            (Ptr{Cvoid},), bonds[i])
        if !hasstereo
            styles[i] = :none
            continue
        end
        iswedge = ccall(
            (:isWedge, libcoordgen), Bool, (Ptr{Cvoid},), bonds[i])
        isrev = ccall(
            (:isReversed, libcoordgen), Bool, (Ptr{Cvoid},), bonds[i])
        if iswedge
            styles[i] = isrev ? :revup : :up
        else
            styles[i] = isrev ? :revdown : :down
        end
    end
    # TODO: keep wave bond in SDFile
    return coords, styles
end

function coordgen!(mol::MolGraph)
    coords, styles = coordgen(mol)
    set_state!(mol, :v_coords2d, coords)
    set_state!(mol, :e_single_bond_style, styles)
    return
end