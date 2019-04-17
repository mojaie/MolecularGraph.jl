#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    sdfilewriter


function atomblock(io::IO, mol::GraphMol)
    for (i, atom) in enumerate(nodeattrs(mol))
        # TODO: SMILES coords
        (x, y, z) = atom.coords
        xyzsym = @sprintf "%10.4f%10.4f%10.4f %-3s" x y z string(atom.symbol)
        println(io, "$(xyzsym) 0  0  0  0  0  0  0  0  0  0  0  0")
    end
    return
end


function bondblock(io::IO, mol::GraphMol)
    for (i, (u, v)) in enumerate(edgesiter(mol))
        bond = edgeattr(mol, i)
        # TODO: SmilesBond
        uv = @sprintf "%3d%3d%3d%3d  0  0  0" u v bond.order bond.notation
        println(io, uv)
    end
    return
end


function propertyblock(io::IO, mol::GraphMol)
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


function datablock(io::IO, mol::GraphMol)
    for (key, val) in mol.attributes
        println(io, "> <$(string(key))>")
        println(io, string(val))
        println(io, "")
    end
    println(io, raw"$$$$")
    return
end


function molblock(io::IO, mol::GraphMol)
    # TODO: 2D coords for SmilesMol
    println(io)
    println(io, "MolecularGraph.jl version $(version())")
    println(io)
    ncnt = nodecount(mol)
    ecnt = edgecount(mol)
    header = @sprintf "%3d%3d  0  0  0  0  0  0  0  0999 V2000" ncnt ecnt
    println(io, header)
    atomblock(io, mol)
    bondblock(io, mol)
    propertyblock(io, mol)
    println(io, "M  END")
    datablock(io, mol)
    return
end


"""
    sdfilewriter(io::IO, mols)
    sdfilewriter(filename::AbstractString, mols)

Write molecule data to the output stream as a SDFile format file.
"""
sdfilewriter(io::IO, mols) = molblock.((io,), mols)
sdfilewriter(filename::AbstractString, mols
    ) = sdfilewriter(open(path, "w"), mols)
