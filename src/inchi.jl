#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

mutable struct inchi_Output
    szInChI::Cstring
    szAuxInfo::Cstring
    szMessage::Cstring
    szLog::Cstring
    inchi_Output() = new(C_NULL, C_NULL, C_NULL, C_NULL)
end

mutable struct inchi_InputINCHI
    # the caller is responsible for the data allocation and deallocation
    szInChI::Cstring       # InChI ASCIIZ string to be converted to a strucure
    szOptions::Cstring     # InChI options: space-delimited; each is preceded by
                           # '/' or '-' depending on OS and compiler
end

mutable struct inchi_Input
    # the caller is responsible for the data allocation and deallocation
    atom::Ptr{Cvoid}            # array of num_atoms elements
    stereo0D::Ptr{Cvoid}        # array of num_stereo0D 0D stereo elements or NULL
    szOptions::Cstring          # InChI options: space-delimited; each is preceded by
                                # '/' or '-' depending on OS and compiler
    num_atoms::Cshort           # number of atoms in the structure < MAX_ATOMS
    num_stereo0D::Cshort        # number of 0D stereo elements
end

mutable struct inchi_InputEx
    # the caller is responsible for the data allocation and deallocation
    atom::Ptr{Cvoid}            # array of num_atoms elements
    stereo0D::Ptr{Cvoid}        # array of num_stereo0D 0D stereo elements or NULL
    szOptions::Cstring          # InChI options: space-delimited; each is preceded by
                                # '/' or '-' depending on OS and compiler
    num_atoms::Cshort           # number of atoms in the structure < MAX_ATOMS
    num_stereo0D::Cshort        # number of 0D stereo elements
    polymer::Ptr{Cvoid}
    v3000::Ptr{Cvoid}
end

mutable struct inchi_OutputStruct
    # the caller is responsible for the data allocation and deallocation
    atom::Ptr{Cvoid}        # array of num_atoms elements
    stereo0D::Ptr{Cvoid}    # array of num_stereo0D 0D stereo elements or NULL
    num_atoms::Cshort;      # number of atoms in the structure < MAX_ATOMS
    num_stereo0D::Cshort    # number of 0D stereo elements
    szMessage::Cstring      # Error/warning ASCIIZ message
    szLog::Cstring          # log-file ASCIIZ string, contains a human-readable list
                            # of recognized options and possibly an Error/warn message
    warningflags::NTuple{2, NTuple{2, Clong}} # warnings, see INCHIDIFF in inchicmp.h                        */
                            # [x][y]:
                            # x=0 => Reconnected if present in InChI
                            #         otherwise Disconnected/Normal
                            # x=1 => Disconnected layer if Reconnected layer is present
                            # y=1 => Main layer or Mobile-H
                            # y=0 => Fixed-H layer

    inchi_OutputStruct() = new(C_NULL, C_NULL, 0, 0, C_NULL, C_NULL, ((0, 0), (0, 0)))
end

mutable struct inchi_OutputStructEx
    # the caller is responsible for the data allocation and deallocation
    atom::Ptr{Cvoid}        # array of num_atoms elements
    stereo0D::Ptr{Cvoid}    # array of num_stereo0D 0D stereo elements or NULL
    num_atoms::Cshort;      # number of atoms in the structure < MAX_ATOMS
    num_stereo0D::Cshort    # number of 0D stereo elements
    szMessage::Cstring      # Error/warning ASCIIZ message
    szLog::Cstring          # log-file ASCIIZ string, contains a human-readable list
                            # of recognized options and possibly an Error/warn message
    warningflags::NTuple{2, NTuple{2, Clong}} # warnings, see INCHIDIFF in inchicmp.h                        */
                            # [x][y]:
                            # x=0 => Reconnected if present in InChI
                            #         otherwise Disconnected/Normal
                            # x=1 => Disconnected layer if Reconnected layer is present
                            # y=1 => Main layer or Mobile-H
                            # y=0 => Fixed-H layer
    polymer::Ptr{Cvoid}
    v3000::Ptr{Cvoid}

    inchi_OutputStructEx() = new(C_NULL, C_NULL, 0, 0, C_NULL, C_NULL, ((0, 0), (0, 0)), C_NULL, C_NULL)
end

function opt_array(options::String)
    isempty(options) ? SubString{String}[] : lstrip.(split(options, ' ', keepempty = false), Ref(['-', '/']))
end

function opt_string(options::Vector{<:AbstractString})
    join((Sys.iswindows() ? '/' : '-') .* options, ' ')
end

function unsafe_info(cs::Cstring, title::String = "")
    if cs != C_NULL
        @info string(title, isempty(title) ? "" : ": ", unsafe_string(cs))
    end
end

function report_output(output, verbose::Bool)
    if output.szInChI == C_NULL || verbose
        output.szInChI == C_NULL && @info "InChI error with \$inchi or \$options"
        unsafe_info(output.szMessage, "message")
        unsafe_info(output.szMessage, "log")
    end
end

"""
    inchi(molblock::String; options::String = "", verbose::Bool = false) -> Union{String,Nothing}
    inchi(mol::MolGraph; options::String = "", verbose::Bool = false) -> Union{String,Nothing}

Generate InChI string from molblock string or molecule.

Options, e.g. "SNon" for 'no stereo information' are specified in https://github.com/mojaie/libinchi/blob/master/INCHIBASE/src/inchiapi.h
"""
function inchi(molblock::String; options::String = "", verbose::Bool = false)
    # support the correct options format depending on OS
    # add a timeout of 60s per molecule
    opts = opt_array(options)
    any(occursin.(r"^Wm?\d+$", opts)) || push!(opts, "W60")
    options = opt_string(opts)

    output = inchi_Output()
    @ccall libinchi.MakeINCHIFromMolfileText(
        molblock::Cstring, options::Cstring, output::Ref{inchi_Output})::Int32
    report_output(output, verbose)

    res = output.szInChI == C_NULL ? nothing : unsafe_string(output.szInChI)

    # Free string buffers allocated by MakeINCHIFromMolfileText
    @ccall libinchi.FreeINCHI(output::Ref{inchi_Output})::Cvoid
    return res
end

inchi(mol::SimpleMolGraph; options::String = "", verbose = false) = inchi(printv2mol(mol); options, verbose)


"""
    inchikey(inchi::String) -> Union{String,Nothing}
    inchikey(mol::MolGraph) -> Union{String,Nothing}

Generate InChI key from InChI string or molecule
"""
function inchikey(inchi::Union{String,Nothing})
    inchi === nothing && return nothing
    ikeybuf = pointer(Vector{UInt8}(undef, 256))
    # TODO: need extra buffer?
    xtra1buf = pointer(Vector{UInt8}(undef, 256))
    xtra2buf = pointer(Vector{UInt8}(undef, 256))
    @ccall libinchi.GetINCHIKeyFromINCHI(
        inchi::Cstring, 1::Int32, 1::Int32,
        ikeybuf::Cstring, xtra1buf::Cstring, xtra2buf::Cstring)::Int32
    return unsafe_string(ikeybuf)
end

inchikey(mol::SimpleMolGraph) = inchikey(inchi(mol))

"""
    inchitosdf(inchi::String; options::String = "") -> Union{String,Nothing}

Generate sdf string from inchi string, `options` are specified in https://github.com/mojaie/libinchi/blob/master/INCHI_BASE/src/inchi_api.h
"""
function inchitosdf(inchi::String; options::String = "", verbose::Bool = false)
    # support the correct options format depending on OS
    opts = opt_array(options)
    # add a timeout of 60s per molecule
    any(occursin.(r"^Wm?\d+$", opts)) || push!(opts, "W60")
    # switch output to sdf format
    "OutputSDF" ∈ opts || push!(opts, "OutputSDF")
    options = opt_string(opts)

    structure = inchi_OutputStructEx()
    inchi_input = inchi_InputINCHI(Base.unsafe_convert(Cstring, inchi), Base.unsafe_convert(Cstring, options))
    @ccall libinchi.GetStructFromINCHIEx(
        inchi_input::Ref{inchi_InputINCHI}, structure::Ref{inchi_OutputStructEx})::Int32

    input = inchi_InputEx(
        structure.atom,
        structure.stereo0D,
        Base.unsafe_convert(Cstring, options),
        structure.num_atoms,
        structure.num_stereo0D,
        structure.polymer,
        structure.v3000
    )
    output = inchi_Output()
    @ccall libinchi.GetINCHIEx(
        input::Ref{inchi_InputEx}, output::Ref{inchi_Output})::Cint
    report_output(output, verbose)

    res = output.szInChI == C_NULL ? nothing : unsafe_string(output.szInChI)

    # Free buffers allocated by GetStructFromINCHIEx and GetINCHI
    @ccall libinchi.FreeStructFromINCHIEx(structure::Ref{inchi_OutputStructEx})::Cvoid
    @ccall libinchi.FreeINCHI(output::Ref{inchi_Output})::Cvoid

    return res
end

"""
    function inchitomol(inchi::String; options = "", verbose = false)

Generate molecule from inchi string, `options` are specified in https://github.com/mojaie/libinchi/blob/master/INCHI_BASE/src/inchi_api.h
"""
function inchitomol(inchi::String; options = "", verbose = false)
    inchitosdf(inchi; options, verbose) |> IOBuffer |> sdftomol
end