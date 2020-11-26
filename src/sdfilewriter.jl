#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    printv2mol, printv2sdf, sdfilewriter


function printv2atoms(io::IO, mol::SDFile)
    for (i, atom) in enumerate(nodeattrs(mol))
        x, y, z = atom.coords
        xyzsym = @sprintf "%10.4f%10.4f%10.4f %-3s" x y z string(atom.symbol)
        println(io, "$(xyzsym) 0  0  0  0  0  0  0  0  0  0  0  0")
    end
    return
end

function printv2atoms(io::IO, mol::SMILES, coords)
    for (i, atom) in enumerate(nodeattrs(mol))
        x, y = coords[i, 1:2]
        z = 0.0
        xyzsym = @sprintf "%10.4f%10.4f%10.4f %-3s" x y z string(atom.symbol)
        println(io, "$(xyzsym) 0  0  0  0  0  0  0  0  0  0  0  0")
    end
    return
end


function printv2bonds(io::IO, mol::SDFile)
    for (i, (u, v)) in enumerate(edgesiter(mol))
        bond = edgeattr(mol, i)
        uv = @sprintf "%3d%3d%3d%3d  0  0  0" u v bond.order bond.notation
        println(io, uv)
    end
    return
end


function printv2bonds(io::IO, mol::SMILES, styles)
    for (i, (u, v)) in enumerate(edgesiter(mol))
        bond = edgeattr(mol, i)
        if styles[i] in (2, 7)
            f, s = (v, u)
            notation = styles[i] - 1
        else
            f, s = (u, v)
            notation = styles[i]
        end
        uv = @sprintf "%3d%3d%3d%3d  0  0  0" f s bond.order notation
        println(io, uv)
    end
    return
end


function printv2properties(io::IO, mol::GraphMol)
    charges = Tuple{Int,Int}[]
    radicals = Tuple{Int,Int}[]
    masses = Tuple{Int,Float64}[]
    for (i, atom) in enumerate(nodeattrs(mol))
        atom.charge == 0 || push!(charges, (i, atom.charge))
        atom.multiplicity == 1 || push!(radicals, (i, atom.multiplicity))
        atom.mass === nothing || push!(masses, (i, atom.mass))
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
    return
end


function printv2data(io::IO, mol::GraphMol)
    for (key, val) in mol.attributes
        println(io, "> <$(string(key))>")
        println(io, string(val))
        println(io, "")
    end
    return
end


function printv2mol(io::IO, mol::GraphMol)
    println(io)
    println(io, "MolecularGraph.jl version $(Util.VERSION)")
    println(io)
    ncnt = nodecount(mol)
    ecnt = edgecount(mol)
    header = @sprintf "%3d%3d  0  0  0  0  0  0  0  0999 V2000" ncnt ecnt
    println(io, header)
    if nodeattrtype(mol) === SmilesAtom
        coords, styles = coordgen(mol)
        printv2atoms(io, mol, coords)
        printv2bonds(io, mol, styles)
    else
        printv2atoms(io, mol)
        printv2bonds(io, mol)
    end
    printv2properties(io, mol)
    println(io, "M  END")
    return
end


function printv2mol(mol::GraphMol)
    buf = IOBuffer(write=true)
    printv2mol(buf, mol)
    res = String(take!(buf))
    close(buf)
    return res
end


function printv2sdf(io::IO, mol::GraphMol)
    printv2mol(io, mol)
    printv2data(io, mol)
    println(io, raw"$$$$")
    return
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
