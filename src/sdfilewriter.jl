#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

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
        io::IO, g::SimpleGraph, atomsymbol::Vector{Symbol}, coords::Vector{Point2d})
    for i in vertices(g)
        x, y = coords[i][1:2]
        z = 0.0  # TODO: keep 3D
        xyzsym = @sprintf "%10.4f%10.4f%10.4f %-3s" x y z string(atomsymbol[i])
        println(io, "$(xyzsym) 0  0  0  0  0  0  0  0  0  0  0  0")
    end
end


function printv2atoms(
        io::IO, g::SimpleGraph, atomsymbol::Vector{Symbol}, coords::Vector{Point3d})
    for i in vertices(g)
        x, y, z = coords[i][1:3]
        xyzsym = @sprintf "%10.4f%10.4f%10.4f %-3s" x y z string(atomsymbol[i])
        println(io, "$(xyzsym) 0  0  0  0  0  0  0  0  0  0  0  0")
    end
end


function printv2bonds(
        io::IO, g::SimpleGraph, bondorder::Vector{Int}, styles::Vector{Symbol})  # 2D
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


printv2mol(io::IO, mol::SimpleMolGraph; kwargs...
    ) = printv2mol(io::IO, mol::SimpleMolGraph, vproptype(mol), eproptype(mol); kwargs...)

function printv2mol(
        io::IO, mol::SimpleMolGraph, V::Type{<:StandardAtom}, E::Type{<:StandardBond}
        ; givebackhydrogen=true)
    # stereospecific hydrogens for aesthetics of fused rings
    # TODO: this may add unnecessary hydrogens and unexpectedly call coordgen!
    # TODO: HydrogenatedAtom should be implemented to keep hydrogen coordinates
    if givebackhydrogen && has_prop(mol, :stereocenter) && !isempty(mol[:stereocenter])
        ringcount = ring_count(mol)
        imph = implicit_hydrogens(mol)
        mol_ = copy(mol)
        for center in keys(mol_.gprops.stereocenter)
            if imph[center] == 1 && ringcount[center] > 1
                add_vertex!(mol_, V(symbol=:H))
                add_edge!(mol_, center, nv(mol_), E())
            end
        end
        mol = mol_
    end
    # write
    program = "MGjlv$(string(VERSION.major)[end])$(string(VERSION.minor)[end-1:end])"
    datetime = Dates.format(Dates.now(), "mmddyyHHMM")
    println(io)
    println(io, "  $(program)$(datetime)2D            ")
    println(io)
    ncnt = nv(mol)
    ecnt = ne(mol)
    header = @sprintf "%3d%3d  0  0  0  0  0  0  0  0999 V2000" ncnt ecnt
    println(io, header)
    bondorder = bond_order(mol)  # dispatch update
    if has_coords3d(mol)
        printv2atoms(io, mol.graph, atom_symbol(mol), coords3d(mol))
        printv2bonds(io, mol.graph, bondorder)
    elseif has_coords2d(mol)
        printv2atoms(io, mol.graph, atom_symbol(mol), coords2d(mol))
        printv2bonds(io, mol.graph, bondorder, draw2d_bond_style(mol))
    else  # Generate coords
        coords, styles = coordgen(mol)
        printv2atoms(io, mol.graph, atom_symbol(mol), coords)
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
function sdfilewriter(io::IO, mols; writer=printv2sdf)
    cnt = length(writer.((io,), mols))
    @info "$(cnt) records exported."
end

function sdfilewriter(filename::AbstractString, mols; kwargs...)
    open(filename, "w") do io
        sdfilewriter(io, mols; kwargs...)
    end
end
