#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

"""
    LibMolGraphJL

C-callable shared library of MolecularGraph.jl built with JuliaC (`juliac --output-lib --trim=safe`).

Molecules are passed as JSON text (`JSON.json(mol)`). Options are passed as JSON objects
(an empty string or `{}` uses the defaults). Strings returned by the library are allocated
by `malloc` and should be released by `free_string`.

On error, functions return a sentinel value (`-1` for integers, `NaN` for floats, an empty string
or an empty molecule for strings) and the error message can be obtained by `last_error`.
"""
module LibMolGraphJL

import Cairo  # activates MolecularGraph's CairoExt (drawpng)
import JSON
import MolecularGraph
using MolecularGraph:
    SMILESMolGraph, SDFMolGraph, SMARTSMolGraph, QueryMolGraph, mol_from_json,
    extract_largest_component!, remove_all_hydrogens!, printv2mol, printv2sdf,
    tdmcis, tdmces, tdmcis_constraints, tdmces_constraints,
    maximum_common_subgraph, maximum_clique_mg
using MolecularGraph.Graphs: SimpleGraph, Edge, nv, ne
import StructUtils


# Options (concretely typed, parsed by `JSON.parse(json, T)`)

StructUtils.@kwarg struct MolOptions
    extract_largest_component::Bool = false
    remove_all_hydrogens::Bool = false
end

StructUtils.@kwarg struct DrawOptions
    viewbox::Bool = true  # SVG only
    show_carbon::String = "simple"  # "simple", "terminal" or "all"
    bgcolor::String = "#FFF"
    bgopacity::Float64 = 1.0
end

StructUtils.@kwarg struct MatchOptions
    stereo::Bool = false
    timeout::Int = 0  # seconds, 0: no timeout
end

StructUtils.@kwarg struct MCSOptions
    diameter::Int = 0  # 0: no distance cutoff
    tolerance::Int = 0
    timeout::Int = 10  # seconds
end

function parse_options(::Type{T}, options::Cstring) where T
    str = unsafe_string(options)
    isempty(strip(str)) && return T()
    return JSON.parse(str, T)
end


# Strings and errors

# Base64 encoding (Base64.base64encode is not compilable with `juliac --trim`)
const BASE64_CHARS = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"

function base64encode(data::Vector{UInt8})
    out = UInt8[]
    sizehint!(out, 4 * cld(length(data), 3))
    for i in 1:3:length(data)
        n = min(3, length(data) - i + 1)
        b = UInt32(data[i]) << 16
        n > 1 && (b |= UInt32(data[i + 1]) << 8)
        n > 2 && (b |= UInt32(data[i + 2]))
        push!(out, BASE64_CHARS[(b >> 18) & 0x3f + 1], BASE64_CHARS[(b >> 12) & 0x3f + 1])
        push!(out, n > 1 ? BASE64_CHARS[(b >> 6) & 0x3f + 1] : UInt8('='))
        push!(out, n > 2 ? BASE64_CHARS[b & 0x3f + 1] : UInt8('='))
    end
    return String(out)
end

const LAST_ERROR = Ref{String}("")

function julia_strdup(s::String)::Cstring
    # Allocate String on the C side
    n = sizeof(s) + 1
    p = Base.Libc.malloc(n)
    p == C_NULL && error("malloc failed")
    unsafe_copyto!(Ptr{UInt8}(p), pointer(s), n)
    return Cstring(p)
end

# Message of an exception without `show` (`show` of arbitrary exceptions is not compilable with trim)
function exception_message(e)
    e isa ErrorException || return string(nameof(typeof(e)))
    m = e.msg
    return m isa String ? m : "error"
end

function set_error!(e)
    LAST_ERROR[] = exception_message(e)
    return
end

"""
    last_error() -> char*

Return the message of the last error (empty if no error has occurred).
"""
Base.@ccallable function last_error()::Cstring
    return julia_strdup(LAST_ERROR[])
end

"""
    free_string(char*)

Release a string returned by the library.
"""
Base.@ccallable function free_string(p::Ptr{UInt8})::Cvoid
    p == C_NULL || Base.Libc.free(p)
    return
end


# Molecule readers (return molecule JSON)

function preprocess!(mol, op::MolOptions)
    op.extract_largest_component && extract_largest_component!(mol)  # default false
    op.remove_all_hydrogens && remove_all_hydrogens!(mol)  # default false
    return mol
end

Base.@ccallable function smilestomol(smiles::Cstring, options::Cstring)::Cstring
    try
        mol = MolecularGraph.smilestomol(unsafe_string(smiles))
        preprocess!(mol, parse_options(MolOptions, options))
        return julia_strdup(JSON.json(mol))
    catch e
        set_error!(e)
        mol = SMILESMolGraph()
        mol[:logs]["error_smiles"] = LAST_ERROR[]
        return julia_strdup(JSON.json(mol))
    end
end

Base.@ccallable function smartstomol(smarts::Cstring)::Cstring
    try
        mol = MolecularGraph.smartstomol(unsafe_string(smarts))
        return julia_strdup(JSON.json(mol))
    catch e
        set_error!(e)
        mol = SMARTSMolGraph()
        mol[:logs]["error_smarts"] = LAST_ERROR[]
        return julia_strdup(JSON.json(mol))
    end
end

Base.@ccallable function sdftomol(sdf::Cstring, options::Cstring)::Cstring
    try
        mol = MolecularGraph.sdftomol(IOBuffer(unsafe_string(sdf)))
        preprocess!(mol, parse_options(MolOptions, options))
        return julia_strdup(JSON.json(mol))
    catch e
        set_error!(e)
        mol = SDFMolGraph()
        mol[:logs]["error_sdfile"] = LAST_ERROR[]
        return julia_strdup(JSON.json(mol))
    end
end


# Properties

Base.@ccallable function vertex_count(mol::Cstring)::Cint
    try
        return Cint(mol_from_json(nv, unsafe_string(mol)))
    catch e
        set_error!(e)
        return Cint(-1)
    end
end

Base.@ccallable function edge_count(mol::Cstring)::Cint
    try
        return Cint(mol_from_json(ne, unsafe_string(mol)))
    catch e
        set_error!(e)
        return Cint(-1)
    end
end

Base.@ccallable function inchikey(mol::Cstring)::Cstring
    try
        ikey = mol_from_json(unsafe_string(mol)) do m
            something(MolecularGraph.inchikey(m), "")
        end
        return julia_strdup(ikey)
    catch e
        set_error!(e)
        return julia_strdup("")
    end
end

Base.@ccallable function standard_weight(mol::Cstring)::Cdouble
    try
        return mol_from_json(m -> MolecularGraph.standard_weight(m, 2), unsafe_string(mol))
    catch e
        set_error!(e)
        return NaN
    end
end


# Writers

Base.@ccallable function molblock(mol::Cstring)::Cstring
    try
        molb = mol_from_json(m -> printv2mol(m; givebackhydrogen=false), unsafe_string(mol))
        return julia_strdup(molb)
    catch e
        set_error!(e)
        return julia_strdup("")
    end
end

Base.@ccallable function sdfmolblock(mol::Cstring)::Cstring
    try
        res = mol_from_json(unsafe_string(mol)) do m
            buf = IOBuffer()
            printv2sdf(buf, m; givebackhydrogen=false)
            String(take!(buf))
        end
        return julia_strdup(res)
    catch e
        set_error!(e)
        return julia_strdup("")
    end
end


# Drawing

Base.@ccallable function drawsvg(mol::Cstring, options::Cstring)::Cstring
    try
        op = parse_options(DrawOptions, options)
        svg = mol_from_json(unsafe_string(mol)) do m
            MolecularGraph.drawsvg(m; viewbox=op.viewbox, show_carbon=Symbol(op.show_carbon),
                bgcolor=op.bgcolor, bgopacity=op.bgopacity)
        end
        return julia_strdup(svg)
    catch e
        set_error!(e)
        return julia_strdup("")
    end
end

"""
    drawpng(mol, width, height, options) -> char*

Return a PNG image encoded in Base64.
"""
Base.@ccallable function drawpng(
        mol::Cstring, width::UInt32, height::UInt32, options::Cstring)::Cstring
    try
        op = parse_options(DrawOptions, options)
        data = mol_from_json(unsafe_string(mol)) do m
            buf = IOBuffer()
            MolecularGraph.drawpng(buf, m, Int(width), Int(height);
                show_carbon=Symbol(op.show_carbon), bgcolor=op.bgcolor, bgopacity=op.bgopacity)
            take!(buf)
        end
        return julia_strdup(base64encode(data))
    catch e
        set_error!(e)
        return julia_strdup("")
    end
end


# Structure match

# A query molecule (SMARTS) cannot be matched against a non-query molecule as `mol1`
is_matchable(m1, m2) = !(m1 isa QueryMolGraph && !(m2 isa QueryMolGraph))

function run_match(f::F, mol1::Cstring, mol2::Cstring, options::Cstring) where F
    op = parse_options(MatchOptions, options)
    mol_from_json(unsafe_string(mol1)) do m1
        mol_from_json(unsafe_string(mol2)) do m2
            is_matchable(m1, m2) || error("query molecule cannot be the first argument")
            op.timeout > 0 ? f(m1, m2; stereo=op.stereo, timeout=op.timeout) : f(m1, m2; stereo=op.stereo)
        end
    end
end

Base.@ccallable function has_exact_match(
        mol1::Cstring, mol2::Cstring, options::Cstring)::Cint
    try
        return Cint(run_match(MolecularGraph.has_exact_match, mol1, mol2, options))
    catch e
        set_error!(e)
        return Cint(-1)
    end
end

Base.@ccallable function has_substruct_match(
        mol1::Cstring, mol2::Cstring, options::Cstring)::Cint
    try
        return Cint(run_match(MolecularGraph.has_substruct_match, mol1, mol2, options))
    catch e
        set_error!(e)
        return Cint(-1)
    end
end


# Maximum common substructure (query molecules are not supported)

function mcs(f::F, mol1::Cstring, mol2::Cstring, options::Cstring) where F
    op = parse_options(MCSOptions, options)
    mol_from_json(unsafe_string(mol1)) do m1
        mol_from_json(unsafe_string(mol2)) do m2
            (m1 isa QueryMolGraph || m2 isa QueryMolGraph) && error("query molecules are not supported")
            if op.diameter > 0
                f(m1, m2; diameter=op.diameter, tolerance=op.tolerance, timeout=op.timeout)
            else
                f(m1, m2; tolerance=op.tolerance, timeout=op.timeout)
            end
        end
    end
end

Base.@ccallable function tdmcis_size(
        mol1::Cstring, mol2::Cstring, options::Cstring)::Cint
    try
        return Cint(mcs((m1, m2; kw...) -> length(tdmcis(m1, m2; kw...)[1]), mol1, mol2, options))
    catch e
        set_error!(e)
        return Cint(-1)
    end
end

Base.@ccallable function tdmces_size(
        mol1::Cstring, mol2::Cstring, options::Cstring)::Cint
    try
        return Cint(mcs((m1, m2; kw...) -> length(tdmces(m1, m2; kw...)[1]), mol1, mol2, options))
    catch e
        set_error!(e)
        return Cint(-1)
    end
end

# Graph-based local similarity (GLS)
function gls(constraints::F, m1, m2; kwargs...) where F
    (nv(m1) == 0 || nv(m2) == 0) && return 0.0
    c1 = constraints(m1; kwargs...)
    c2 = constraints(m2; kwargs...)
    m1max = length(maximum_clique_mg(SimpleGraph(Edge.(c1.pairs)))[1])
    m2max = length(maximum_clique_mg(SimpleGraph(Edge.(c2.pairs)))[1])
    comm = length(maximum_common_subgraph(c1, c2; kwargs...)[1])
    return comm / (m1max + m2max - comm)
end

Base.@ccallable function tdmcis_gls(
        mol1::Cstring, mol2::Cstring, options::Cstring)::Cdouble
    try
        return mcs((m1, m2; kw...) -> gls(tdmcis_constraints, m1, m2; kw...), mol1, mol2, options)
    catch e
        set_error!(e)
        return NaN
    end
end

Base.@ccallable function tdmces_gls(
        mol1::Cstring, mol2::Cstring, options::Cstring)::Cdouble
    try
        return mcs((m1, m2; kw...) -> gls(tdmces_constraints, m1, m2; kw...), mol1, mol2, options)
    catch e
        set_error!(e)
        return NaN
    end
end

end  # module
