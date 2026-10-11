#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module LibMolGraphJL

using Base: unsafe_convert
using Base64: Base64EncodePipe
using Cairo
using Graphs
using JSON
using MolecularGraph


function julia_strdup(s::String)::Cstring
    # Allocate String on the C side
    n = sizeof(s) + 1
    p = Base.Libc.malloc(n)
    p == C_NULL && error("malloc failed")
    unsafe_copyto!(Ptr{UInt8}(p), pointer(s), n)
    return Cstring(p)
end


Base.@ccallable function smilestomol(smiles::Cstring, options::Cstring)::Cstring
    try
        smstr = unsafe_string(smiles)  # explicit copy
        opstr = unsafe_string(options)  # explicit copy
        mol = MolecularGraph.smilestomol(smstr)
        op = JSON.parse(opstr)
        if haskey(op, "extract_largest_component") && op["extract_largest_component"]
            extract_largest_component!(mol)  # default extract_largest_component=false
        end
        if haskey(op, "remove_all_hydrogens") && op["remove_all_hydrogens"]
            remove_all_hydrogens!(mol)  # default remove_all_hydrogens=false
        end
        julia_strdup(JSON.json(mol))
    catch e
        mol = SMILESMolGraph()
        mol[:logs]["error_smiles"] = e.msg
        julia_strdup(JSON.json(mol))
    end
end


Base.@ccallable function smartstomol(smarts::Cstring)::Cstring
    try
        smstr = unsafe_string(smarts)  # explicit copy
        mol = MolecularGraph.smartstomol(smstr)
        julia_strdup(JSON.json(mol))
    catch e
        mol = SMARTSMolGraph()
        mol[:logs]["error_smarts"] = e.msg
        julia_strdup(JSON.json(mol))
    end
end


Base.@ccallable function sdftomol(sdf::Cstring, options::Cstring)::Cstring
    try
        # return empty mol on error
        sdfstr = unsafe_string(sdf)  # explicit copy
        opstr = unsafe_string(options)  # explicit copy
        mol = MolecularGraph.sdftomol(IOBuffer(sdfstr))
        op = JSON.parse(opstr)
        if haskey(op, "extract_largest_component") && op["extract_largest_component"]
            extract_largest_component!(mol)  # default extract_largest_component=false
        end
        if haskey(op, "remove_all_hydrogens") && op["remove_all_hydrogens"]
            remove_all_hydrogens!(mol)  # default remove_all_hydrogens=false
        end
        julia_strdup(JSON.json(mol))
    catch e
        mol = SDFMolGraph()
        mol[:logs]["error_sdfile"] = e.msg
        julia_strdup(JSON.json(mol))
    end
end


Base.@ccallable function vertex_count(mol::Cstring)::Cint
    try
        mstr = unsafe_string(mol)  # explicit copy
        molobj = mol_from_json(mstr)
        nv(molobj)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function edge_count(mol::Cstring)::Cint
    try
        mstr = unsafe_string(mol)  # explicit copy
        molobj = mol_from_json(mstr)
        ne(molobj)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function inchikey(mol::Cstring)::Cstring
    try
        mstr = unsafe_string(mol)  # explicit copy
        molobj = mol_from_json(mstr)
        ikey = MolecularGraph.inchikey(molobj)
        julia_strdup(something(ikey, ""))
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function standard_weight(mol::Cstring)::Cdouble
    try
        mstr = unsafe_string(mol)  # explicit copy
        molobj = mol_from_json(mstr)
        MolecularGraph.standard_weight(molobj, 2)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function molblock(mol::Cstring)::Cstring
    try
        mstr = unsafe_string(mol)  # explicit copy
        molobj = mol_from_json(mstr)
        molb = printv2mol(molobj; givebackhydrogen=false)
        julia_strdup(molb)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function sdfmolblock(mol::Cstring)::Cstring
    try
        mstr = unsafe_string(mol)  # explicit copy
        molobj = mol_from_json(mstr)
        buf = IOBuffer(write=true)
        printv2sdf(buf, molobj; givebackhydrogen=false)
        res = String(take!(buf))
        close(buf)
        julia_strdup(res)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function drawsvg(mol::Cstring, options::Cstring)::Cstring
    try
        mstr = unsafe_string(mol)  # explicit copy
        opstr = unsafe_string(options)  # explicit copy
        molobj = mol_from_json(mstr)
        op = JSON.parse(opstr)
        kwgs = Pair{Symbol,Any}[]
        haskey(op, "viewbox") && push!(kwgs, :viewbox => op["viewbox"])
        haskey(op, "show_carbon") && push!(kwgs, :show_carbon => Symbol(op["show_carbon"]))
        haskey(op, "bgcolor") && push!(kwgs, :bgcolor => op["bgcolor"])
        haskey(op, "bgopacity") && push!(kwgs, :bgopacity => op["bgopacity"])
        svg = MolecularGraph.drawsvg(molobj; kwgs...)
        julia_strdup(svg)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function drawpng(
        mol::Cstring, width::UInt32, height::UInt32, options::Cstring)::Cstring
    try
        mstr = unsafe_string(mol)  # explicit copy
        opstr = unsafe_string(options)  # explicit copy
        molobj = mol_from_json(mstr)
        op = JSON.parse(opstr)
        kwgs = Pair{Symbol,Any}[]
        haskey(op, "show_carbon") && push!(kwgs, :show_carbon => Symbol(op["show_carbon"]))
        haskey(op, "bgcolor") && push!(kwgs, :bgcolor => op["bgcolor"])
        haskey(op, "bgopacity") && push!(kwgs, :bgopacity => op["bgopacity"])
        buf = IOBuffer()
        iob64_encode = Base64EncodePipe(buf)
        MolecularGraph.drawpng(iob64_encode, molobj, Int(width), Int(height); kwgs...)
        close(iob64_encode)
        str = String(take!(buf))
        close(buf)
        julia_strdup(str)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function has_exact_match(
        mol1::Cstring, mol2::Cstring, kwargs::Cstring)::Cint
    try
        mstr1 = unsafe_string(mol1)  # explicit copy
        mstr2 = unsafe_string(mol2)  # explicit copy
        kwstr = unsafe_string(kwargs)  # explicit copy
        molobj1 = mol_from_json(mstr1)
        molobj2 = mol_from_json(mstr2)
        kw = JSON.parse(kwstr)
        MolecularGraph.has_exact_match(molobj1, molobj2; kw...)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function has_substruct_match(
        mol1::Cstring, mol2::Cstring, kwargs::Cstring)::Cint
    try
        mstr1 = unsafe_string(mol1)  # explicit copy
        mstr2 = unsafe_string(mol2)  # explicit copy
        kwstr = unsafe_string(kwargs)  # explicit copy
        molobj1 = mol_from_json(mstr1)
        molobj2 = mol_from_json(mstr2)
        kw = JSON.parse(kwstr)
        MolecularGraph.has_substruct_match(molobj1, molobj2; kw...)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function tdmcis_size(
        mol1::Cstring, mol2::Cstring, kwargs::Cstring)::Cint
    try
        mstr1 = unsafe_string(mol1)  # explicit copy
        mstr2 = unsafe_string(mol2)  # explicit copy
        kwstr = unsafe_string(kwargs)  # explicit copy
        molobj1 = mol_from_json(mstr1)
        molobj2 = mol_from_json(mstr2)
        kw = Dict(Symbol(k) => v for (k, v) in JSON.parse(kwstr))
        length(tdmcis(molobj1, molobj2; kw...)[1])
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function tdmces_size(
        mol1::Cstring, mol2::Cstring, kwargs::Cstring)::Cint
    try
        mstr1 = unsafe_string(mol1)  # explicit copy
        mstr2 = unsafe_string(mol2)  # explicit copy
        kwstr = unsafe_string(kwargs)  # explicit copy
        molobj1 = mol_from_json(mstr1)
        molobj2 = mol_from_json(mstr2)
        kw = Dict(Symbol(k) => v for (k, v) in JSON.parse(kwstr))
        length(tdmces(molobj1, molobj2; kw...)[1])
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function tdmcis_gls(
        mol1::Cstring, mol2::Cstring, kwargs::Cstring)::Cdouble
    try
        mstr1 = unsafe_string(mol1)  # explicit copy
        mstr2 = unsafe_string(mol2)  # explicit copy
        kwstr = unsafe_string(kwargs)  # explicit copy
        molobj1 = mol_from_json(mstr1)
        molobj2 = mol_from_json(mstr2)
        if nv(molobj1) == 0 || nv(molobj2) == 0
            return 0.0
        end
        kw = Dict(Symbol(k) => v for (k, v) in JSON.parse(kwstr))
        mol1_ = tdmcis_constraints(molobj1; kw...)
        mol2_ = tdmcis_constraints(molobj2; kw...)
        m1max = length(maximum_clique(SimpleGraph(Edge.(mol1_.pairs)))[1])
        m2max = length(maximum_clique(SimpleGraph(Edge.(mol2_.pairs)))[1])
        comm = length(maximum_common_subgraph(mol1_, mol2_; kw...)[1])
        comm / (m1max + m2max - comm)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function tdmces_gls(
        mol1::Cstring, mol2::Cstring, kwargs::Cstring)::Cdouble
    try
        mstr1 = unsafe_string(mol1)  # explicit copy
        mstr2 = unsafe_string(mol2)  # explicit copy
        kwstr = unsafe_string(kwargs)  # explicit copy
        molobj1 = mol_from_json(mstr1)
        molobj2 = mol_from_json(mstr2)
        if nv(molobj1) == 0 || nv(molobj2) == 0
            return 0.0
        end
        kwargs = Dict(Symbol(k) => v for (k, v) in JSON.parse(kwstr))
        mol1_ = tdmces_constraints(molobj1; kwargs...)
        mol2_ = tdmces_constraints(molobj2; kwargs...)
        m1max = length(maximum_clique(SimpleGraph(Edge.(mol1_.pairs)))[1])
        m2max = length(maximum_clique(SimpleGraph(Edge.(mol2_.pairs)))[1])
        comm = length(maximum_common_subgraph(mol1_, mol2_; kwargs...)[1])
        comm / (m1max + m2max - comm)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

end  # module