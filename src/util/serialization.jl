
# Note: Unmarshal patch for Set
# Note: Not used. Get rid of sets. Save as vectors
function Unmarshal.unmarshal(::Type{Set{E}},
    parsedJson::Union{Vector, AbstractArray},
    verbose::Bool=false, verboseLvl::Int=0) where E
    if verbose
        prettyPrint(verboseLvl, "Set{\$E}")
        verboseLvl += 1
    end
    Set(Unmarshal.unmarshal(
        E, field, verbose, verboseLvl) for field in parsedJson)
end