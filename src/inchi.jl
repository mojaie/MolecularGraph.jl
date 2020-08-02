#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    inchi, inchikey

using libinchi_jll


mutable struct inchi_Output
    szInChI::Cstring
    szAuxInfo::Cstring
    szMessage::Cstring
    szLog::Cstring
    inchi_Output() = new(C_NULL, C_NULL, C_NULL, C_NULL)
end


"""
    inchi(molblock::String) -> String
    inchi(mol::GraphMol) -> String

Generate InChI string from molblock string or molecule
"""
function inchi(molblock::String)
    output = inchi_Output()
    ccall(
        (:MakeINCHIFromMolfileText, libinchi),
        Int32,
        (Cstring, Cstring, Ref{inchi_Output}),
        molblock, "-W60", output)
    if output.szInChI == C_NULL
        println("InChI error with $(molblock)")
        if output.szMessage != C_NULL
            println("message: ", unsafe_string(output.szMessage))
        end
        if output.szLog != C_NULL
            println("log: ", unsafe_string(output.szLog))
        end
        res = nothing
    else
        res = unsafe_string(output.szInChI)
    end
    # Free string buffers allocated by MakeINCHIFromMolfileText
    ccall(
        (:FreeINCHI, libinchi),
        Cvoid,
        (Ref{inchi_Output},),
        output
    )
    return res
end

inchi(mol::GraphMol) = inchi(printv2mol(mol))


"""
    inchikey(inchi::String) -> String
    inchikey(mol::GraphMol) -> String

Generate InChI key from InChI string or molecule
"""
function inchikey(inchi::Union{String,Nothing})
    inchi === nothing && return nothing
    ikeybuf = pointer(Vector{UInt8}(undef, 256))
    # TODO: need extra buffer?
    xtra1buf = pointer(Vector{UInt8}(undef, 256))
    xtra2buf = pointer(Vector{UInt8}(undef, 256))
    ccall(
        (:GetINCHIKeyFromINCHI, libinchi),
        Int32,
        (Cstring, Int32, Int32, Cstring, Cstring, Cstring),
        inchi, 1, 1, ikeybuf, xtra1buf, xtra2buf)
    return unsafe_string(ikeybuf)
end

inchikey(mol::GraphMol) = inchikey(inchi(mol))
