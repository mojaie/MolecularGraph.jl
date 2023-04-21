#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    has_coords, sdf_coords2d, coords2d, coords3d, coordgen, coordgen!

using coordgenlibs_jll


function has_coords(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_coords2d) && return true
    hasfield(vproptype(mol), :coords) || return false
    for i in vertices(mol)
        isnothing(get_prop(mol, i, :coords)) && return false
    end
    return true
end

function sdf_coords2d(mol::SimpleMolGraph)
    coords = zeros(Float64, nv(mol), 2)
    for i in vertices(mol)
        coords[i, :] = get_prop(mol, i, :coords)[1:2]
    end
    return coords
end

function coords2d(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_coords2d) && return get_cache(mol, :v_coords2d)
    has_coords(mol) || error("no coordinates. use coordgen")
    return sdf_coords2d(mol)
end


function coords3d(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_coords3d) && return get_cache(mol, :v_coords3d)
    has_coords(mol) || error("no coordinates. use coordgen")
    coords = zeros(Float64, nv(mol), 3)
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
function coordgen(g, atomsymbol_, bondorder_, stereocenters, stereobonds)
    minmol = ccall((:getSketcherMinimizer, libcoordgen), Ptr{Cvoid}, ())
    atoms = Ptr{Cvoid}[]
    bonds = Ptr{Cvoid}[]
    er = Dict(e => i for (i, e) in enumerate(edges(g)))

    # Atoms
    for a in atomsymbol_
        atom = ccall(
            (:setAtom, libcoordgen), Ptr{Cvoid},
            (Ptr{Cvoid}, Int), minmol, atomnumber(a)
        )
        push!(atoms, atom)
    end

    # Bonds
    for (i, e) in enumerate(edges(g))
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
            bonds[er[e]], atoms[stereo[1]], atoms[stereo[2]], stereo[3]
        )
    end

    # Optimize
    ccall(
        (:runGenerateCoordinates, libcoordgen), Cvoid,
        (Ptr{Cvoid},), minmol
    )

    # Output
    coords = zeros(Float64, nv(g), 2)
    for i in vertices(g)
        px = ccall(
            (:getAtomX, libcoordgen), Float32, (Ptr{Cvoid},), atoms[i])
        py = ccall(
            (:getAtomY, libcoordgen), Float32, (Ptr{Cvoid},), atoms[i])
        coords[i, :] = [px, py]
    end
    styles = Vector{Symbol}(undef, ne(g))
    for (i, e) in enumerate(edges(g))
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


function coordgen(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    return coordgen(
        mol.graph, atom_symbol(mol), bond_order(mol),
        get_prop(mol, :stereocenter), get_prop(mol, :stereobond)
    )
end


function coordgen!(mol::SimpleMolGraph)
    coords, styles = coordgen(
        mol.graph, atom_symbol(mol), bond_order(mol),
        get_prop(mol, :stereocenter), get_prop(mol, :stereobond)
    )
    set_cache!(mol, :v_coords2d, coords)
    set_cache!(mol, :e_coordgen_bond_style, styles)
end