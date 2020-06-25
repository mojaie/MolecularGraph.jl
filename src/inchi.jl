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

Generate InChI string from molblock
"""
function inchi(molblock::String)
    output = inchi_Output()
    ccall(
        (:MakeINCHIFromMolfileText, libinchi),
        Int32,
        (Cstring, Cstring, Ref{inchi_Output}),
        molblock, "-W60", output)
    res = unsafe_string(output.szInChI)
    # Free string buffers allocated by MakeINCHIFromMolfileText
    ccall(
        (:FreeINCHI, libinchi),
        Cvoid,
        (Ref{inchi_Output},),
        output
    )
    return res
end


"""
    inchi(mol::GraphMol) -> String

Generate InChI string from GraphMol
"""
function inchi(mol::GraphMol)
    buf = IOBuffer(write=true)
    molblock(buf, mol)
    return inchi(String(take!(buf)))
end


"""
    inchikey(inchi::String) -> String

Generate InChI key from InChI
"""
function inchikey(inchi::String)
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
