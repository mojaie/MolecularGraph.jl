#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# Note: just to avoid type piracy like StructUtils.lower(x::Point2d)


function coords_from_sdf!(mol::SimpleMolGraph)
    # Initializer. should be called before stereochem initializers.
    zrange = nv(mol) == 0 ? (0, 0) : extrema(mol[i].coords[3] for i in vertices(mol))
    # Embed 3D coords to 2D
    push!(
        mol[:descriptors].coords2d,
        Coords2d([Point2d(mol[i].coords[1:2]...) for i in vertices(mol)]))
    if zrange[2] - zrange[1] > 0.001  # 3D coords available
        push!(
            mol[:descriptors].coords3d,
            Coords3d([Point3d(mol[i].coords[1:3]...) for i in vertices(mol)]))
    end
    # Bond style in 2D notation (wedges)
    bondorder = [mol[e].order for e in edges(mol)]
    bondnotation = [mol[e].notation for e in edges(mol)]
    isordered = [mol[e].isordered for e in edges(mol)]
    arr = Vector{Symbol}(undef, length(bondorder))
    for i in 1:length(bondorder)
        if bondnotation[i] == 3
            arr[i] = :cis_trans
        elseif bondorder[i] != 1
            arr[i] = :none
        elseif bondnotation[i] == 1
            arr[i] = isordered[i] ? :up : :revup
        elseif bondnotation[i] == 6
            arr[i] = isordered[i] ? :down : :revdown
        elseif bondnotation[i] == 4
            arr[i] = :unspecified
        else
            arr[i] = :none
        end
    end
    push!(mol[:descriptors].draw2d_bond_style, Draw2dBondStyle(arr))
end


function coords2d(mol::SimpleMolGraph, i::Integer)
    if !has_descriptor(mol, :coords2d)
        error("2D coordinates not available in this molecule type.")
    end
    cds = get_descriptor(mol, :coords2d)
    if length(cds) < i
        error("No available coords found. Run coordgen!(mol) to generate coords")
    elseif nv(mol) != length(cds[i])
        error("Coords size mismatch. Run coordgen!(mol) to regenerate coords")
    end
    return cds[i]
end
coords2d(mol::SimpleMolGraph) = coords2d(mol, 1)

function has_coords2d(mol::SimpleMolGraph)
    has_descriptor(mol, :coords2d) || return false
    isempty(get_descriptor(mol, :coords2d)) && return false
    length(get_descriptor(mol, :coords2d)[1]) == nv(mol) || return false
    return true
end


function coords3d(mol::SimpleMolGraph, i::Integer)
    if !has_descriptor(mol, :coords3d)
        error("3D coordinates not available in this molecule type.")
    end
    cds = get_descriptor(mol, :coords3d)
    if length(cds) < i
        error("No available 3D coords found.")
    elseif nv(mol) != length(cds[i])
        error("Coords size mismatch.")
    end
    return cds[i]
end
coords3d(mol::SimpleMolGraph) = coords3d(mol, 1)

function has_coords3d(mol::SimpleMolGraph)
    has_descriptor(mol, :coords3d) || return false
    isempty(get_descriptor(mol, :coords3d)) && return false
    length(get_descriptor(mol, :coords3d)[1]) == nv(mol) || return false
    return true
end


draw2d_bond_style(mol::SimpleMolGraph, i::Integer
    ) = get_descriptor(mol, :draw2d_bond_style)[i]
draw2d_bond_style(mol::SimpleMolGraph) = draw2d_bond_style(mol, 1)


function update_coords!(mol::SimpleMolGraph)
    # TODO: just remap nodes if still all existing vertices have coords.
    empty!(mol[:descriptors].coords2d)
    empty!(mol[:descriptors].coords3d)
    empty!(mol[:descriptors].draw2d_bond_style)
    if hasfield(vproptype(mol), :coords)
        if !any(isnothing(mol[i].coords) for i in vertices(mol))
            coords_from_sdf!(mol)
        end
    else
        coordgen!(mol)
    end
end



"""
    coordgen(mol::MolGraph) -> Tuple{Array{Int,1},Array{Int,1}}

Generate 2D coordinates by using Schrodinger's coordgenlibs.

This will returns a tuple of `coords` and `styles` arrays. `coords` is a size(n, 2) matrix
where n is atom count, which stores 2D coordinates (x, y) of each atoms.
`styles` is a size e vector of wedge notation of stereobond, where e is bond count.
"""
coordgen(mol::SimpleMolGraph) = coordgen(
    mol.graph, atom_number(mol), bond_order(mol),
    mol[:stereocenter], mol[:stereobond]
)

function coordgen(
        g::SimpleGraph{T}, atomnum::Vector{Int}, bondorder_::Vector{Int},
        stereocenters::StereocenterMap{T}, stereobonds::StereobondMap{T}) where T
    minmol = @ccall libcoordgen.getSketcherMinimizer()::Ptr{Cvoid}
    atoms = Ptr{Cvoid}[]
    bonds = Ptr{Cvoid}[]
    ernk = edge_rank(g)

    # Atoms
    for a in atomnum
        a = a > 0 ? a : 9  # F as a placeholder of virtual atom
        atom = @ccall libcoordgen.setAtom(
            minmol::Ptr{Cvoid}, a::Cint)::Ptr{Cvoid}
        push!(atoms, atom)
    end

    # Bonds
    for (i, e) in enumerate(edges(g))
        bond = @ccall libcoordgen.setBond(
            minmol::Ptr{Cvoid}, atoms[src(e)]::Ptr{Cvoid},
            atoms[dst(e)]::Ptr{Cvoid}, bondorder_[i]::Cint)::Ptr{Cvoid}
        push!(bonds, bond)
    end

    @ccall libcoordgen.assignBondsAndNeighbors(minmol::Ptr{Cvoid})::Cvoid

    # Stereocenter
    for (n, c) in stereocenters
        @ccall libcoordgen.setStereoCenter(
            atoms[n]::Ptr{Cvoid}, atoms[c.lookingFrom]::Ptr{Cvoid},
            atoms[c.first]::Ptr{Cvoid}, atoms[c.second]::Ptr{Cvoid}, c.isclockwise::Cint)::Cvoid
    end

    # Stereobond
    for (e, b) in stereobonds
        @ccall libcoordgen.setStereoBond(
            bonds[ernk[e]]::Ptr{Cvoid}, atoms[b.first]::Ptr{Cvoid},
            atoms[b.second]::Ptr{Cvoid}, b.is_cis::Cint
        )::Cvoid
    end

    # Optimize
    @ccall libcoordgen.runGenerateCoordinates(minmol::Ptr{Cvoid})::Cvoid

    # Output
    coords = Vector{Point2d}(undef, nv(g))
    for i in vertices(g)
        px = @ccall libcoordgen.getAtomX(atoms[i]::Ptr{Cvoid})::Cfloat
        py = @ccall libcoordgen.getAtomY(atoms[i]::Ptr{Cvoid})::Cfloat
        coords[i] = Point2d(px, py)
    end
    coords = coords * 0.0165  # scale to practical coords scale (bond length = 0.825)
    styles = Vector{Symbol}(undef, ne(g))
    for (i, e) in enumerate(edges(g))
        hasstereo = @ccall libcoordgen.hasStereochemistryDisplay(
            bonds[i]::Ptr{Cvoid})::Bool
        if !hasstereo
            styles[i] = :none
            continue
        end
        iswedge = @ccall libcoordgen.isWedge(bonds[i]::Ptr{Cvoid})::Bool
        isrev = @ccall libcoordgen.isReversed(bonds[i]::Ptr{Cvoid})::Bool
        if iswedge
            styles[i] = isrev ? :revup : :up
        else
            styles[i] = isrev ? :revdown : :down
        end
    end
    # TODO: keep wave bond in SDFile
    return Coords2d(coords), Draw2dBondStyle(styles)
end

function coordgen!(mol::SimpleMolGraph)
    # Initialize
    empty!(mol[:descriptors].coords2d)
    empty!(mol[:descriptors].draw2d_bond_style)
    # TODO: unspecified stereochem in SMILES
    coords, styles = coordgen(
        mol.graph, atom_number(mol), bond_order(mol),
        mol[:stereocenter], mol[:stereobond]
    )
    push!(mol[:descriptors].coords2d, coords)
    push!(mol[:descriptors].draw2d_bond_style, styles)
end
