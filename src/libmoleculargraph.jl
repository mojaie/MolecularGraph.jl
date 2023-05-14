#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

Base.@ccallable function smilestomol(smiles::Ptr{UInt8})::Ptr{UInt8}
    return try
        mol = smilestomol(unsafe_string(smiles))
        buf = IOBuffer(write=true)
        JSON.print(buf, to_dict(mol))
        return pointer(buf.data)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function smartstomol(smarts::Ptr{UInt8})::Ptr{UInt8}
    return try
        mol = smartstomol(unsafe_string(smarts))
        buf = IOBuffer(write=true)
        JSON.print(buf, to_dict(mol))
        return pointer(buf.data)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function sdftomol(sdf::Ptr{UInt8})::Ptr{UInt8}
    return try
        mol = sdftomol(IOBuffer(unsafe_string(sdf)))
        buf = IOBuffer(write=true)
        JSON.print(buf, to_dict(mol))
        return pointer(buf.data)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function inchikey(mol::Ptr{UInt8})::Ptr{UInt8}
    return try
        mol = MolGraph(JSON.parse(unsafe_string(mol)))
        ikey = inchikey(mol)
        buf = IOBuffer(write=true)
        print(buf, ikey)
        return pointer(buf.data)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function standard_weight(mol::Ptr{UInt8})::Cdouble
    return try
        mol = MolGraph(JSON.parse(unsafe_string(mol)))
        return standard_weight(mol, 2)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function has_exact_match(
        mol1::Ptr{UInt8}, mol2::Ptr{UInt8}, kwargs::Ptr{UInt8})::Cint
    return try
        mol1 = MolGraph(JSON.parse(unsafe_string(mol1)))
        mol2 = MolGraph(JSON.parse(unsafe_string(mol2)))
        kwargs = JSON.parse(unsafe_string(kwargs))
        return has_exact_match(mol1, mol2; kwargs...)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function has_substruct_match(
        mol1::Ptr{UInt8}, mol2::Ptr{UInt8}, kwargs::Ptr{UInt8})::Cint
    return try
        mol1 = MolGraph(JSON.parse(unsafe_string(mol1)))
        mol2 = MolGraph(JSON.parse(unsafe_string(mol2)))
        kwargs = JSON.parse(unsafe_string(kwargs))
        return has_substruct_match(mol1, mol2; kwargs...)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function tcmcis(
        mol1::Ptr{UInt8}, mol2::Ptr{UInt8}, kwargs::Ptr{UInt8})::Cint
    return try
        mol1 = MolGraph(JSON.parse(unsafe_string(mol1)))
        mol2 = MolGraph(JSON.parse(unsafe_string(mol2)))
        kwargs = Dict(Symbol(k) => v for (k, v) in JSON.parse(unsafe_string(kwargs)))
        return size(tcmcis(mol1, mol2; kwargs...))
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function tcmces(
        mol1::Ptr{UInt8}, mol2::Ptr{UInt8}, kwargs::Ptr{UInt8})::Cint
    return try
        mol1 = MolGraph(JSON.parse(unsafe_string(mol1)))
        mol2 = MolGraph(JSON.parse(unsafe_string(mol2)))
        kwargs = Dict(Symbol(k) => v for (k, v) in JSON.parse(unsafe_string(kwargs)))
        return size(tcmces(mol1, mol2; kwargs...))
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end