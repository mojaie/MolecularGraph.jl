#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

using JSON
export
    smilestomol, sdftomol,
    inchikey, standardweight,
    hasexactmatch, hassubstructmatch, tcmcis, tcmces


Base.@ccallable function smilestomol(smiles::Ptr{UInt8})::Ptr{UInt8}
    return try
        mol = smilestomol(unsafe_string(smiles))
        buf = IOBuffer(write=true)
        JSON.print(buf, todict(mol))
        return pointer(buf.data)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function sdftomol(sdf::Ptr{UInt8})::Ptr{UInt8}
    return try
        mol = sdftomol(split(unsafe_string(sdf), "\n"))
        buf = IOBuffer(write=true)
        JSON.print(buf, todict(mol))
        return pointer(buf.data)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function inchikey(mol::Ptr{UInt8})::Ptr{UInt8}
    return try
        mol = graphmol(JSON.parse(unsafe_string(mol)))
        ikey = inchikey(mol)
        buf = IOBuffer(write=true)
        print(buf, ikey)
        return pointer(buf.data)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function standardweight(mol::Ptr{UInt8})::Cdouble
    return try
        mol = graphmol(JSON.parse(unsafe_string(mol)))
        return standardweight(Float64, mol)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function hasexactmatch(
        mol1::Ptr{UInt8}, mol2::Ptr{UInt8}, kwargs::Ptr{UInt8})::Cint
    return try
        mol1 = graphmol(JSON.parse(unsafe_string(mol1)))
        mol2 = graphmol(JSON.parse(unsafe_string(mol2)))
        kwargs = JSON.parse(unsafe_string(kwargs))
        return hasexactmatch(mol1, mol2; kwargs...)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function hassubstructmatch(
        mol1::Ptr{UInt8}, mol2::Ptr{UInt8}, kwargs::Ptr{UInt8})::Cint
    return try
        mol1 = graphmol(JSON.parse(unsafe_string(mol1)))
        mol2 = graphmol(JSON.parse(unsafe_string(mol2)))
        kwargs = JSON.parse(unsafe_string(kwargs))
        return hassubstructmatch(mol1, mol2; kwargs...)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function tcmcis(
        mol1::Ptr{UInt8}, mol2::Ptr{UInt8}, kwargs::Ptr{UInt8})::Cint
    return try
        mol1 = graphmol(JSON.parse(unsafe_string(mol1)))
        mol2 = graphmol(JSON.parse(unsafe_string(mol2)))
        kwargs = Dict(Symbol(k) => v for (k, v) in JSON.parse(unsafe_string(kwargs)))
        return size(tcmcis(mol1, mol2; kwargs...))
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function tcmces(
        mol1::Ptr{UInt8}, mol2::Ptr{UInt8}, kwargs::Ptr{UInt8})::Cint
    return try
        mol1 = graphmol(JSON.parse(unsafe_string(mol1)))
        mol2 = graphmol(JSON.parse(unsafe_string(mol2)))
        kwargs = Dict(Symbol(k) => v for (k, v) in JSON.parse(unsafe_string(kwargs)))
        return size(tcmces(mol1, mol2; kwargs...))
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end