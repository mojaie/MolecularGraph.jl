#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    AnnotationError,
    OperationError

import Base: showerror


struct AnnotationError <: Exception
    msg::String
end

showerror(io::IO, e::AnnotationError) = print(io, e.msg)


struct OperationError <: Exception
    msg::String
end

showerror(io::IO, e::OperationError) = print(io, e.msg)
