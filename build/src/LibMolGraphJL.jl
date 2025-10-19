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


Base.@ccallable function smilestomol(smiles::Cstring, options::Cstring)::Cstring
    try
        mol = MolecularGraph.smilestomol(unsafe_string(smiles))
        op = JSON.parse(unsafe_string(options))
        if haskey(op, "extract_largest_component") && op["extract_largest_component"]
            extract_largest_component!(mol)  # default extract_largest_component=false
        end
        if haskey(op, "remove_all_hydrogens") && op["remove_all_hydrogens"]
            remove_all_hydrogens!(mol)  # default remove_all_hydrogens=false
        end
        unsafe_convert(Cstring, JSON.json(mol))
    catch e
        mol = SMILESMolGraph()
        mol[:logs]["error_smiles"] = e.msg
        unsafe_convert(Cstring, JSON.json(mol))
    end
end


Base.@ccallable function smartstomol(smarts::Cstring)::Cstring
    try
        mol = MolecularGraph.smartstomol(unsafe_string(smarts))
        unsafe_convert(Cstring, JSON.json(mol))
    catch e
        mol = SMARTSMolGraph()
        mol[:logs]["error_smarts"] = e.msg
        unsafe_convert(Cstring, JSON.json(mol))
    end
end


Base.@ccallable function sdftomol(sdf::Cstring, options::Cstring)::Cstring
    try
        # return empty mol on error
        mol = MolecularGraph.sdftomol(IOBuffer(unsafe_string(sdf)))
        op = JSON.parse(unsafe_string(options))
        if haskey(op, "extract_largest_component") && op["extract_largest_component"]
            extract_largest_component!(mol)  # default extract_largest_component=false
        end
        if haskey(op, "remove_all_hydrogens") && op["remove_all_hydrogens"]
            remove_all_hydrogens!(mol)  # default remove_all_hydrogens=false
        end
        unsafe_convert(Cstring, JSON.json(mol))
    catch e
        mol = SDFMolGraph()
        mol[:logs]["error_sdfile"] = e.msg
        unsafe_convert(Cstring, JSON.json(mol))
    end
end


Base.@ccallable function vertex_count(mol::Cstring)::Cint
    try
        molobj = mol_from_json(unsafe_string(mol))
        nv(molobj)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function edge_count(mol::Cstring)::Cint
    try
        molobj = mol_from_json(unsafe_string(mol))
        ne(molobj)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function inchikey(mol::Cstring)::Cstring
    try
        molobj = mol_from_json(unsafe_string(mol))
        ikey = MolecularGraph.inchikey(molobj)
        unsafe_convert(Cstring, something(ikey, ""))
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function standard_weight(mol::Cstring)::Cdouble
    try
        molobj = mol_from_json(unsafe_string(mol))
        MolecularGraph.standard_weight(molobj, 2)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function molblock(mol::Cstring)::Cstring
    try
        molobj = mol_from_json(unsafe_string(mol))
        unsafe_convert(Cstring, printv2mol(molobj; givebackhydrogen=false))
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function sdfmolblock(mol::Cstring)::Cstring
    try
        molobj = mol_from_json(unsafe_string(mol))
        buf = IOBuffer(write=true)
        printv2sdf(buf, molobj; givebackhydrogen=false)
        res = String(take!(buf))
        close(buf)
        unsafe_convert(Cstring, res)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function drawsvg(mol::Cstring, options::Cstring)::Cstring
    try
        molobj = mol_from_json(unsafe_string(mol))
        op = JSON.parse(unsafe_string(options))
        kwgs = Pair{Symbol,Any}[]
        haskey(op, "viewbox") && push!(kwgs, :viewbox => op["viewbox"])
        haskey(op, "show_carbon") && push!(kwgs, :show_carbon => Symbol(op["show_carbon"]))
        haskey(op, "bgcolor") && push!(kwgs, :bgcolor => op["bgcolor"])
        haskey(op, "bgopacity") && push!(kwgs, :bgopacity => op["bgopacity"])
        svg = MolecularGraph.drawsvg(molobj; kwgs...)
        unsafe_convert(Cstring, svg)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function drawpng(
        mol::Cstring, width::UInt32, height::UInt32, options::Cstring)::Cstring
    try
        molobj = mol_from_json(unsafe_string(mol))
        op = JSON.parse(unsafe_string(options))
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
        unsafe_convert(Cstring, str)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function has_exact_match(
        mol1::Cstring, mol2::Cstring, kwargs::Cstring)::Cint
    try
        mol1 = mol_from_json(unsafe_string(mol1))
        mol2 = mol_from_json(unsafe_string(mol2))
        kwargs = JSON.parse(unsafe_string(kwargs))
        MolecularGraph.has_exact_match(mol1, mol2; kwargs...)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function has_substruct_match(
        mol1::Cstring, mol2::Cstring, kwargs::Cstring)::Cint
    try
        mol1 = mol_from_json(unsafe_string(mol1))
        mol2 = mol_from_json(unsafe_string(mol2))
        kwargs = JSON.parse(unsafe_string(kwargs))
        MolecularGraph.has_substruct_match(mol1, mol2; kwargs...)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function tdmcis_size(
        mol1::Cstring, mol2::Cstring, kwargs::Cstring)::Cint
    try
        mol1 = mol_from_json(unsafe_string(mol1))
        mol2 = mol_from_json(unsafe_string(mol2))
        kwargs = Dict(Symbol(k) => v for (k, v) in JSON.parse(unsafe_string(kwargs)))
        length(tdmcis(mol1, mol2; kwargs...)[1])
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function tdmces_size(
        mol1::Cstring, mol2::Cstring, kwargs::Cstring)::Cint
    try
        mol1 = mol_from_json(unsafe_string(mol1))
        mol2 = mol_from_json(unsafe_string(mol2))
        kwargs = Dict(Symbol(k) => v for (k, v) in JSON.parse(unsafe_string(kwargs)))
        length(tdmces(mol1, mol2; kwargs...)[1])
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function tdmcis_gls(
        mol1::Cstring, mol2::Cstring, kwargs::Cstring)::Cdouble
    try
        mol1 = mol_from_json(unsafe_string(mol1))
        mol2 = mol_from_json(unsafe_string(mol2))
        if nv(mol1) == 0 || nv(mol2) == 0
            return 0.0
        end
        kwargs = Dict(Symbol(k) => v for (k, v) in JSON.parse(unsafe_string(kwargs)))
        mol1_ = tdmcis_constraints(mol1; kwargs...)
        mol2_ = tdmcis_constraints(mol2; kwargs...)
        m1max = length(maximum_clique(SimpleGraph(Edge.(mol1_.pairs)))[1])
        m2max = length(maximum_clique(SimpleGraph(Edge.(mol2_.pairs)))[1])
        comm = length(maximum_common_subgraph(mol1_, mol2_; kwargs...)[1])
        comm / (m1max + m2max - comm)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end


Base.@ccallable function tdmces_gls(
        mol1::Cstring, mol2::Cstring, kwargs::Cstring)::Cdouble
    try
        mol1 = mol_from_json(unsafe_string(mol1))
        mol2 = mol_from_json(unsafe_string(mol2))
        if ne(mol1) == 0 || ne(mol2) == 0
            return 0.0
        end
        kwargs = Dict(Symbol(k) => v for (k, v) in JSON.parse(unsafe_string(kwargs)))
        mol1_ = tdmces_constraints(mol1; kwargs...)
        mol2_ = tdmces_constraints(mol2; kwargs...)
        m1max = length(maximum_clique(SimpleGraph(Edge.(mol1_.pairs)))[1])
        m2max = length(maximum_clique(SimpleGraph(Edge.(mol2_.pairs)))[1])
        comm = length(maximum_common_subgraph(mol1_, mol2_; kwargs...)[1])
        comm / (m1max + m2max - comm)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

end  # module