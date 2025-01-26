#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

using Base: unsafe_convert
using Base64

export
    vertex_count, edge_count,
    molblock, sdfmolblock,
    tdmcis_size, tdmces_size, tdmcis_gls, tdmces_gls


Base.@ccallable function smilestomol(smiles::Cstring, options::Cstring)::Cstring
    try
        mol = smilestomol(unsafe_string(smiles))
        op = JSON.parse(unsafe_string(options))
        if haskey(op, "extract_largest_component") && op["extract_largest_component"]
            extract_largest_component!(mol)  # default extract_largest_component=false
        end
        if haskey(op, "remove_all_hydrogens") && op["remove_all_hydrogens"]
            remove_all_hydrogens!(mol)  # default remove_all_hydrogens=false
        end
        return unsafe_convert(Cstring, JSON.json(to_dict(mol)))
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function smartstomol(smarts::Cstring)::Cstring
    return try
        mol = smartstomol(unsafe_string(smarts))
        return unsafe_convert(Cstring, JSON.json(to_dict(mol)))
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function sdftomol(sdf::Cstring, options::Cstring)::Cstring
    return try
        # return empty mol on error
        mol = iterate(sdfilereader(IOBuffer(unsafe_string(sdf)); unsupported=:ignore))[1]
        op = JSON.parse(unsafe_string(options))
        if haskey(op, "extract_largest_component") && op["extract_largest_component"]
            extract_largest_component!(mol)  # default extract_largest_component=false
        end
        if haskey(op, "remove_all_hydrogens") && op["remove_all_hydrogens"]
            remove_all_hydrogens!(mol)  # default remove_all_hydrogens=false
        end
        return unsafe_convert(Cstring, JSON.json(to_dict(mol)))
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function vertex_count(mol::Cstring)::Cint
    return try
        molobj = MolGraph(JSON.parse(unsafe_string(mol)))
        return nv(molobj)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function edge_count(mol::Cstring)::Cint
    return try
        molobj = MolGraph(JSON.parse(unsafe_string(mol)))
        return ne(molobj)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function inchikey(mol::Cstring)::Cstring
    return try
        molobj = MolGraph(JSON.parse(unsafe_string(mol)))
        ikey = inchikey(molobj)
        return unsafe_convert(Cstring, something(ikey, ""))
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function standard_weight(mol::Cstring)::Cdouble
    return try
        molobj = MolGraph(JSON.parse(unsafe_string(mol)))
        return standard_weight(molobj, 2)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function molblock(mol::Cstring)::Cstring
    return try
        molobj = MolGraph(JSON.parse(unsafe_string(mol)))
        return unsafe_convert(Cstring, printv2mol(molobj))
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function sdfmolblock(mol::Cstring)::Cstring
    return try
        molobj = MolGraph(JSON.parse(unsafe_string(mol)))
        buf = IOBuffer(write=true)
        printv2sdf(buf, molobj)
        res = String(take!(buf))
        close(buf)
        return unsafe_convert(Cstring, res)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function drawsvg(mol::Cstring)::Cstring
    return try
        molobj = MolGraph(JSON.parse(unsafe_string(mol)))
        svg = drawsvg(molobj)
        return unsafe_convert(Cstring, svg)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function drawpng(
        mol::Cstring, width::UInt32, height::UInt32)::Cstring
    return try
        molobj = MolGraph(JSON.parse(unsafe_string(mol)))
        buf = IOBuffer()
        iob64_encode = Base64EncodePipe(buf)
        drawpng(iob64_encode, molobj, Int(width), Int(height))
        close(iob64_encode)
        str = String(take!(buf))
        close(buf)
        return unsafe_convert(Cstring, str)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function has_exact_match(
        mol1::Cstring, mol2::Cstring, kwargs::Cstring)::Cint
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
        mol1::Cstring, mol2::Cstring, kwargs::Cstring)::Cint
    return try
        mol1 = MolGraph(JSON.parse(unsafe_string(mol1)))
        mol2 = MolGraph(JSON.parse(unsafe_string(mol2)))
        kwargs = JSON.parse(unsafe_string(kwargs))
        return has_substruct_match(mol1, mol2; kwargs...)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function tdmcis_size(
        mol1::Cstring, mol2::Cstring, kwargs::Cstring)::Cint
    return try
        mol1 = MolGraph(JSON.parse(unsafe_string(mol1)))
        mol2 = MolGraph(JSON.parse(unsafe_string(mol2)))
        kwargs = Dict(Symbol(k) => v for (k, v) in JSON.parse(unsafe_string(kwargs)))
        return length(tdmcis(mol1, mol2; kwargs...)[1])
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function tdmces_size(
        mol1::Cstring, mol2::Cstring, kwargs::Cstring)::Cint
    return try
        mol1 = MolGraph(JSON.parse(unsafe_string(mol1)))
        mol2 = MolGraph(JSON.parse(unsafe_string(mol2)))
        kwargs = Dict(Symbol(k) => v for (k, v) in JSON.parse(unsafe_string(kwargs)))
        return length(tdmces(mol1, mol2; kwargs...)[1])
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function tdmcis_gls(
        mol1::Cstring, mol2::Cstring, kwargs::Cstring)::Cdouble
    return try
        mol1 = MolGraph(JSON.parse(unsafe_string(mol1)))
        mol2 = MolGraph(JSON.parse(unsafe_string(mol2)))
        if nv(mol1) == 0 || nv(mol2) == 0
            return 0.0
        end
        kwargs = Dict(Symbol(k) => v for (k, v) in JSON.parse(unsafe_string(kwargs)))
        mol1_ = tdmcis_constraints(mol1; kwargs...)
        mol2_ = tdmcis_constraints(mol2; kwargs...)
        m1max = length(maximum_clique(SimpleGraph(Edge.(mol1_.pairs)))[1])
        m2max = length(maximum_clique(SimpleGraph(Edge.(mol2_.pairs)))[1])
        comm = length(maximum_common_subgraph(mol1_, mol2_; kwargs...)[1])
        return comm / (m1max + m2max - comm)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function tdmces_gls(
        mol1::Cstring, mol2::Cstring, kwargs::Cstring)::Cdouble
    return try
        mol1 = MolGraph(JSON.parse(unsafe_string(mol1)))
        mol2 = MolGraph(JSON.parse(unsafe_string(mol2)))
        if ne(mol1) == 0 || ne(mol2) == 0
            return 0.0
        end
        kwargs = Dict(Symbol(k) => v for (k, v) in JSON.parse(unsafe_string(kwargs)))
        mol1_ = tdmces_constraints(mol1; kwargs...)
        mol2_ = tdmces_constraints(mol2; kwargs...)
        m1max = length(maximum_clique(SimpleGraph(Edge.(mol1_.pairs)))[1])
        m2max = length(maximum_clique(SimpleGraph(Edge.(mol2_.pairs)))[1])
        comm = length(maximum_common_subgraph(mol1_, mol2_; kwargs...)[1])
        return comm / (m1max + m2max - comm)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end
