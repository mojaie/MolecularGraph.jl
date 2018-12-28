#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    GraphMolError,
    MolParseError,
    AnnotationError,
    OperationError

import Base: showerror


abstract type GraphMolError <: Exception end


struct MolParseError <: GraphMolError
    msg::String
end


struct AnnotationError <: GraphMolError
    msg::String
end


struct OperationError <: GraphMolError
    msg::String
end

showerror(io::IO, e::GraphMolError) = show(io, e.msg)
