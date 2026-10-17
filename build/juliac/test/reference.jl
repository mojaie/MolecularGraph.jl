#
# Compare the results of libmolgraphjl (written by python_test.py) with those of the
# LibMolGraphJL module running on Julia.
#
#   julia --project=build/juliac build/juliac/test/reference.jl results.json
#

using LibMolGraphJL
import JSON

function call(func::String, args::Vector)
    f = getfield(LibMolGraphJL, Symbol(func))
    cargs = Any[a isa String ? a : UInt32(a) for a in args]
    ptrs = Any[a isa String ? Base.unsafe_convert(Cstring, a) : a for a in cargs]
    res = GC.@preserve cargs f(ptrs...)
    res isa Cstring || return res
    s = unsafe_string(res)
    LibMolGraphJL.free_string(Ptr{UInt8}(res))
    return s
end

# Mask the timestamp in the header of molfiles
normalize(x) = x isa String ? replace(x, r"(MGjlv\d{3})\d{10}" => s"\1**********") : x

function main(path)
    data = JSON.parse(read(path, String))
    failed = 0
    for case in data["cases"]
        name = case["name"]
        expected = data["results"][name]
        actual = call(case["func"], collect(case["args"]))
        actual isa AbstractFloat && isnan(actual) && (actual = nothing)
        if isequal(normalize(actual), normalize(expected))
            println("ok: ", name)
        else
            failed += 1
            println("NG: ", name, "\n  lib:   ", first(repr(expected), 200),
                "\n  julia: ", first(repr(actual), 200))
        end
    end
    println(failed == 0 ? "all $(length(data["cases"])) cases matched" : "$failed cases differ")
    return failed == 0 ? 0 : 1
end

exit(main(ARGS[1]))
