#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Annotation, AnnotationArray


abstract type Annotation end

struct AnnotationArray
    data::Matrix
    keys::Vector{Symbol}
end

Base.iterate(
    arr::AnnotationArray, state=1
) = state > length(arr) ? nothing : (arr.data[state, :], state + 1)

Base.IteratorEltype(::Type{AnnotationArray}) = Base.EltypeUnknown()
Base.length(arr::AnnotationArray) = size(arr.data, 1)
