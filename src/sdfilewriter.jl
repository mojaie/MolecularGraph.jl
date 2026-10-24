#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# Parsed once at load time (a format string parsed at runtime is not compilable with `juliac --trim`)
const SDF_DATEFORMAT = Dates.DateFormat("mmddyyHHMM")

const BOND_STYLE_TO_SDF = Dict(
    :none => 0,
    :up => 1,
    :revup => 1,
    :down => 6,
    :revdown => 6,
    :unspecified => 4,
    :cis_trans => 3
)


function printv2atoms(
        io::IO, g::SimpleGraph, atomsymbol::Vector{Symbol}, coords::Coords2d, valences::Vector{Int}=Int[])
    for i in vertices(g)
        x, y = coords[i][1:2]
        z = 0.0  # TODO: keep 3D
        xyzsym = @sprintf "%10.4f%10.4f%10.4f %-3s" x y z string(atomsymbol[i])
        val = isempty(valences) ? 0 : valences[i]
        @printf(io, "%s%2d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d\n", xyzsym, 0, 0, 0, 0, 0, val, 0, 0, 0, 0, 0, 0)
    end
end


function printv2atoms(
        io::IO, g::SimpleGraph, atomsymbol::Vector{Symbol}, coords::Coords3d, valences::Vector{Int}=Int[])
    for i in vertices(g)
        x, y, z = coords[i][1:3]
        xyzsym = @sprintf "%10.4f%10.4f%10.4f %-3s" x y z string(atomsymbol[i])
        val = isempty(valences) ? 0 : valences[i]
        @printf(io, "%s%2d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d\n", xyzsym, 0, 0, 0, 0, 0, val, 0, 0, 0, 0, 0, 0)
    end
end


function printv2bonds(
        io::IO, g::SimpleGraph, bondorder::Vector{Int}, styles::Draw2dBondStyle)  # 2D
    for (i, e) in enumerate(edges(g))
        u, v = styles[i] in (:revup, :revdown) ? (dst(e), src(e)) : (src(e), dst(e))
        uv = @sprintf "%3d%3d%3d%3d  0  0  0" u v bondorder[i] BOND_STYLE_TO_SDF[styles[i]]
        println(io, uv)
    end
end

function printv2bonds(io::IO, g::SimpleGraph, bondorder::Vector{Int})  # 3D
    for (i, e) in enumerate(edges(g))
        uv = @sprintf "%3d%3d%3d%3d  0  0  0" src(e) dst(e) bondorder[i] 0
        println(io, uv)
    end
end


function printv2properties(io::IO, mol::SimpleMolGraph)
    charges = Tuple{Int,Int}[]
    radicals = Tuple{Int,Int}[]
    masses = Tuple{Int,Float64}[]
    for i in vertices(mol)
        atom_charge(mol[i]) == 0 || push!(charges, (i, atom_charge(mol[i])))
        multiplicity(mol[i]) == 1 || push!(radicals, (i, multiplicity(mol[i])))
        isotope(mol[i]) == 0 || push!(masses, (i, isotope(mol[i])))
    end
    if !isempty(charges)
        head = @sprintf "M  CHG%3d" length(charges)
        print(io, head)
        for (i, chg) in charges
            rcd = @sprintf " %3d %3d" i chg
            print(io, rcd)
        end
        println(io)
    end
    if !isempty(radicals)
        head = @sprintf "M  RAD%3d" length(radicals)
        print(io, head)
        for (i, rad) in radicals
            rcd = @sprintf " %3d %3d" i rad
            print(io, rcd)
        end
        println(io)
    end
    if !isempty(masses)
        head = @sprintf "M  ISO%3d" length(masses)
        print(io, head)
        for (i, iso) in masses
            rcd = @sprintf " %3d %3d" i iso
            print(io, rcd)
        end
        println(io)
    end
end


function printv2data(io::IO, mol::SimpleMolGraph)
    for (key, val) in mol.gprops.metadata
        println(io, "> <$(string(key))>")
        println(io, string(val))
        println(io, "")
    end
end


# Return a copy of the molecule with explicit hydrogens at stereocenters in fused rings.
# If 2D coordinates are available, the hydrogens are placed at the largest angular gap of
# the neighbors and the wedges of the stereocenters are put on the hydrogens.
function giveback_stereo_hydrogens(mol::SimpleMolGraph, V::Type, E::Type)
    ringcount = ring_count(mol)
    imph = implicit_hydrogens(mol)
    centers = [c for c in keys(mol[:stereocenter]) if imph[c] == 1 && ringcount[c] > 1]
    isempty(centers) && return mol
    has2d = has_coords2d(mol)
    mol_ = copy(mol)
    desc = mol_.gprops.descriptors
    for center in sort(centers)
        if has2d
            for coords in desc.coords2d
                push!(coords.coords, hydrogen_position(coords, center, neighbors(mol_, center)))
            end
        end
        add_vertex!(mol_, V(symbol=:H))
        add_edge!(mol_, center, nv(mol_), E())
        # wedges are re-assigned in the update (preferably on the hydrogen)
        has2d && discard_center_wedges!(mol_, center)
    end
    return mol_
end

function hydrogen_position(coords::Coords2d, center::Integer, nbrs)
    c = coords[center]
    angs = sort([atan(coords[n][2] - c[2], coords[n][1] - c[1]) for n in nbrs])
    gaps = [angs[2:end] .- angs[1:end-1]; angs[1] + 2pi - angs[end]]
    i = argmax(gaps)
    dir = angs[i] + gaps[i] / 2
    len = sum(sqrt((coords[n][1] - c[1])^2 + (coords[n][2] - c[2])^2) for n in nbrs) / length(nbrs)
    return Point2d(c[1] + len * cos(dir), c[2] + len * sin(dir))
end

function discard_center_wedges!(mol::SimpleMolGraph, center::Integer)
    ernk = edge_rank(mol.graph)
    for styles in mol.gprops.descriptors.draw2d_bond_style
        length(styles) == ne(mol) || continue
        for n in neighbors(mol, center)
            i = edge_rank(ernk, center, n)
            st = styles[i]
            if (center < n && st in (:up, :down)) || (center > n && st in (:revup, :revdown))
                styles[i] = :none
            end
        end
    end
    return
end


"""
    printv2mol(io::IO, mol::SimpleMolGraph; givebackhydrogen=true)
    printv2mol(mol::SimpleMolGraph; givebackhydrogen=true) -> String

Write the molecule as a molfile (CTAB V2000) block to the output stream, or return it as a string.

If `givebackhydrogen` is true, hydrogens at stereocenters in fused rings are added for
drawing (except for molecules with 3D coordinates). If 2D coordinates are available, they
are kept and the hydrogens are placed with wedges; otherwise 2D coordinates are generated by
[`coordgen`](@ref).
Available for molecules with `StandardAtom` and `StandardBond`.
"""
printv2mol(io::IO, mol::SimpleMolGraph; kwargs...
    ) = printv2mol(io::IO, mol::SimpleMolGraph, vproptype(mol), eproptype(mol); kwargs...)

function valence_flag(intrinsic_valence::Int, valence::Int)
    intrinsic_valence < 0 && return 0
    return valence == 0 ? 15 : valence > 14 ? 0 : valence
end

function printv2mol(
        io::IO, mol::SimpleMolGraph, V::Type{<:StandardAtom}, E::Type{<:StandardBond}
        ; givebackhydrogen=true)
    # stereospecific hydrogens for aesthetics of fused rings (not for 3D coordinates)
    if (givebackhydrogen && !has_coords3d(mol)
            && haskey(mol, :stereocenter) && !isempty(mol[:stereocenter]))
        mol = giveback_stereo_hydrogens(mol, V, E)
    end

    # Valence column of the atom line for atoms with a fixed number of hydrogens (0: not set, 15: zero valence)
    valences = if hasfield(V, :intrinsic_valence)
        vals = valence(mol)
        Int[valence_flag(mol.vprops[i].intrinsic_valence, vals[i]) for i in vertices(mol)]
    else
        Int[]
    end

    # write
    # 8 characters: the last digit of the major version and the last 2 digits of the minor version
    # (zero-padded without lpad, which is not compilable with juliac --trim)
    minor = string(VERSION.minor)
    minor2 = length(minor) == 1 ? "0" * minor : minor[end-1:end]
    program = "MGjlv$(string(VERSION.major)[end])$(minor2)"
    datetime = Dates.format(Dates.now(), SDF_DATEFORMAT)
    dimension = has_coords3d(mol) ? "3D" : "2D"

    println(io)
    println(io, "  $(program)$(datetime)$(dimension)            ")
    println(io)
    ncnt = nv(mol)
    ecnt = ne(mol)
    header = @sprintf "%3d%3d  0  0  0  0  0  0  0  0999 V2000" ncnt ecnt
    println(io, header)
    bondorder = bond_order(mol)  # dispatch update
    if has_coords3d(mol)
        printv2atoms(io, mol.graph, atom_symbol(mol), coords3d(mol), valences)
        printv2bonds(io, mol.graph, bondorder)
    elseif has_coords2d(mol)
        printv2atoms(io, mol.graph, atom_symbol(mol), coords2d(mol), valences)
        printv2bonds(io, mol.graph, bondorder, draw2d_bond_style(mol))
    else  # Generate coords
        coords, styles = coordgen(mol)
        printv2atoms(io, mol.graph, atom_symbol(mol), coords, valences)
        printv2bonds(io, mol.graph, bondorder, styles)
    end
    printv2properties(io, mol)
    println(io, "M  END")
end


function printv2mol(mol::SimpleMolGraph; kwargs...)
    buf = IOBuffer(write=true)
    printv2mol(buf, mol; kwargs...)
    res = String(take!(buf))
    close(buf)
    return res
end


"""
    printv2sdf(io::IO, mol::SimpleMolGraph; givebackhydrogen=true)

Write the molecule as an SDFile record (molfile block, data fields and `\$\$\$\$`) to the output stream.
"""
function printv2sdf(io::IO, mol::SimpleMolGraph; kwargs...)
    printv2mol(io, mol; kwargs...)
    printv2data(io, mol)
    println(io, raw"$$$$")
end


"""
    sdfilewriter(io::IO, mols)
    sdfilewriter(filename::AbstractString, mols)

Write molecule data to the output stream as a SDFile format file.
"""
function sdfilewriter(io::IO, mols; writer=printv2sdf, kwargs...)
    cnt = length(writer.((io,), mols; kwargs...))
    @info "$(cnt) records exported."
end

function sdfilewriter(filename::AbstractString, mols; kwargs...)
    open(filename, "w") do io
        sdfilewriter(io, mols; kwargs...)
    end
end
