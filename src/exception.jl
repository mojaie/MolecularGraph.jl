#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    MolecularGraphError,
    MolParseError,
    AnnotationError,
    OperationError

import Base: showerror


abstract type MolecularGraphError <: Exception end


struct MolParseError <: MolecularGraphError
    msg::String
end


struct AnnotationError <: MolecularGraphError
    msg::String
end


struct OperationError <: MolecularGraphError
    msg::String
end

showerror(io::IO, e::MolecularGraphError) = show(io, e.msg)
