#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    vertex_count, edge_count,
    tdmcis_size, tdmces_size, tdmcis_gls, tdmces_gls,
    # deprecated
    smilestosvg, sdftosvg,  
    tdmcis_tanimoto, tdmces_tanimoto, tdmcis_dist, tdmces_dist, tdmces_gls_batch


Base.@ccallable function smilestomol(smiles::Ptr{UInt8}, options::Ptr{UInt8})::Ptr{UInt8}
    try
        mol = smilestomol(unsafe_string(smiles))
        op = JSON.parse(unsafe_string(options))
        if haskey(op, "extract_largest_component") && op["extract_largest_component"]
            extract_largest_component!(mol)  # default extract_largest_component=false
        end
        if haskey(op, "remove_all_hydrogens") && op["remove_all_hydrogens"]
            remove_all_hydrogens!(mol)  # default remove_all_hydrogens=false
        end
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


Base.@ccallable function sdftomol(sdf::Ptr{UInt8}, options::Ptr{UInt8})::Ptr{UInt8}
    return try
        # return empty mol on error
        mol = iterate(sdfilereader(IOBuffer(unsafe_string(sdf)); unsupported=:log))[1]
        op = JSON.parse(unsafe_string(options))
        if haskey(op, "extract_largest_component") && op["extract_largest_component"]
            extract_largest_component!(mol)  # default extract_largest_component=false
        end
        if haskey(op, "remove_all_hydrogens") && op["remove_all_hydrogens"]
            remove_all_hydrogens!(mol)  # default remove_all_hydrogens=false
        end
        buf = IOBuffer(write=true)
        JSON.print(buf, to_dict(mol))
        return pointer(buf.data)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function vertex_count(mol::Ptr{UInt8})::Cint
    return try
        mol = MolGraph(JSON.parse(unsafe_string(mol)))
        return nv(mol)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function edge_count(mol::Ptr{UInt8})::Cint
    return try
        mol = MolGraph(JSON.parse(unsafe_string(mol)))
        return ne(mol)
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


Base.@ccallable function drawsvg(mol::Ptr{UInt8})::Ptr{UInt8}
    return try
        mol = MolGraph(JSON.parse(unsafe_string(mol)))
        svg = drawsvg(mol)
        buf = IOBuffer(write=true)
        print(buf, svg)
        return pointer(buf.data)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function drawpng(res::Ptr{UInt8}, mol::Ptr{UInt8}, width::UInt32, height::UInt32)::Cint
    return try
        mol = MolGraph(JSON.parse(unsafe_string(mol)))
        buf = IOBuffer(write=true)
        drawpng(buf, mol, Int(width), Int(height))
        for i = 1:buf.size
            unsafe_store!(res, buf.data[i], i)
        end
        return buf.size
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


Base.@ccallable function tdmcis_size(
        mol1::Ptr{UInt8}, mol2::Ptr{UInt8}, kwargs::Ptr{UInt8})::Cint
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
        mol1::Ptr{UInt8}, mol2::Ptr{UInt8}, kwargs::Ptr{UInt8})::Cint
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
        mol1::Ptr{UInt8}, mol2::Ptr{UInt8}, kwargs::Ptr{UInt8})::Cdouble
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
        mol1::Ptr{UInt8}, mol2::Ptr{UInt8}, kwargs::Ptr{UInt8})::Cdouble
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



# deprecated

Base.@ccallable function tdmces_gls_batch(query::Ptr{UInt8})::Ptr{UInt8}
    # query = [[ids1], [ids2], [queries1], [queries2], {option}, int]
    return try
        ids1, ids2, queries1, queries2, option, thld = JSON.parse(unsafe_string(query))
        kwargs = Dict(Symbol(k) => v for (k, v) in option)
        res = Tuple{Int,Int,Float64}[]
        for (u, v, q1, q2) in zip(ids1, ids2, queries1, queries2)
            mol1 = MolGraph(q1)
            mol2 = MolGraph(q2)
            if ne(mol1) == 0 || ne(mol2) == 0
                continue
            end
            mol1_ = tdmces_constraints(mol1; kwargs...)
            mol2_ = tdmces_constraints(mol2; kwargs...)
            m1max = length(maximum_clique(SimpleGraph(Edge.(mol1_.pairs)))[1])
            m2max = length(maximum_clique(SimpleGraph(Edge.(mol2_.pairs)))[1])
            comm = length(maximum_common_subgraph(mol1_, mol2_; kwargs...)[1])
            score = comm / (m1max + m2max - comm)
            if score >= thld
                push!(res, (u, v, score))
            end
        end
        buf = IOBuffer(write=true)
        JSON.print(buf, res)
        return pointer(buf.data)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function smilestosvg(smiles::Ptr{UInt8})::Ptr{UInt8}
    return try
        mol = smilestomol(unsafe_string(smiles))
        svg = drawsvg(mol)
        buf = IOBuffer(write=true)
        print(buf, svg)
        return pointer(buf.data)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function sdftosvg(sdf::Ptr{UInt8})::Ptr{UInt8}
    return try
        mol = sdftomol(IOBuffer(unsafe_string(sdf)))
        svg = drawsvg(mol)
        buf = IOBuffer(write=true)
        print(buf, svg)
        return pointer(buf.data)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function tdmcis_tanimoto(
        mol1::Ptr{UInt8}, mol2::Ptr{UInt8}, kwargs::Ptr{UInt8})::Cdouble
    return try
        mol1 = MolGraph(JSON.parse(unsafe_string(mol1)))
        mol2 = MolGraph(JSON.parse(unsafe_string(mol2)))
        if nv(mol1) == 0 || nv(mol2) == 0
            return 0.0
        end
        kwargs = Dict(Symbol(k) => v for (k, v) in JSON.parse(unsafe_string(kwargs)))
        comm = length(tdmcis(mol1, mol2; kwargs...)[1])
        return comm / (nv(mol1) + nv(mol2) - comm)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function tdmces_tanimoto(
        mol1::Ptr{UInt8}, mol2::Ptr{UInt8}, kwargs::Ptr{UInt8})::Cdouble
    return try
        mol1 = MolGraph(JSON.parse(unsafe_string(mol1)))
        mol2 = MolGraph(JSON.parse(unsafe_string(mol2)))
        if ne(mol1) == 0 || ne(mol2) == 0
            return 0.0
        end
        kwargs = Dict(Symbol(k) => v for (k, v) in JSON.parse(unsafe_string(kwargs)))
        comm = length(tdmces(mol1, mol2; kwargs...)[1])
        return comm / (ne(mol1) + ne(mol2) - comm)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function tdmcis_dist(
        mol1::Ptr{UInt8}, mol2::Ptr{UInt8}, kwargs::Ptr{UInt8})::Cint
    return try
        mol1 = MolGraph(JSON.parse(unsafe_string(mol1)))
        mol2 = MolGraph(JSON.parse(unsafe_string(mol2)))
        kwargs = Dict(Symbol(k) => v for (k, v) in JSON.parse(unsafe_string(kwargs)))
        comm = length(tdmcis(mol1, mol2; kwargs...)[1])
        return nv(mol1) + nv(mol2) - comm * 2
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function tdmces_dist(
        mol1::Ptr{UInt8}, mol2::Ptr{UInt8}, kwargs::Ptr{UInt8})::Cint
    return try
        mol1 = MolGraph(JSON.parse(unsafe_string(mol1)))
        mol2 = MolGraph(JSON.parse(unsafe_string(mol2)))
        kwargs = Dict(Symbol(k) => v for (k, v) in JSON.parse(unsafe_string(kwargs)))
        comm = length(tdmcis(mol1, mol2; kwargs...)[1])
        return ne(mol1) + ne(mol2) - comm * 2
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end