# import Pkg
# Pkg.add("YAML")
import YAML

const CONVENTIONAL_WEIGHT = Dict{Int,Float64}(
    1 => 1.008,
    3 => 6.94,
    5 => 10.81,
    6 => 12.011,
    7 => 14.007,
    8 => 15.999,
    12 => 24.305,
    14 => 28.085,
    16 => 32.06,
    17 => 35.45,
    35 => 79.904,
    81 => 204.38
)

const path = joinpath(
    dirname(@__FILE__), "../assets/raw/nist_atomic_weights20200410.txt")
const dest = joinpath(
    dirname(@__FILE__), "../assets/const/atomicweights.yaml")
const revdest = joinpath(
    dirname(@__FILE__), "../assets/const/symboltonumber.yaml")

function run()
    data = []
    revdict = Dict{String,Int}()
    lines = eachline(path)
    state = nothing
    next = iterate(lines, state)
    while true
        rcd = Dict()
        iso = Dict()
        while true
            (line, state) = next
            line == "" && break
            k, v = split(line, " = ")
            if k == "Atomic Number"
                rcd["Number"] = parse(Int, v)
            elseif k == "Atomic Symbol"
                rcd["Symbol"] = v
            elseif k == "Mass Number"
                iso["Number"] = parse(Int, v)
            elseif k == "Relative Atomic Mass"
                m = match(r"(.*?)\((.*?)#?\)", v)
                value = parse(Float64, m[1])
                unc = parse(Float64, string(
                    replace(m[1], r"[0-9]" => "0")[1:end-length(m[2])], m[2]))
                iso["Mass"] = value
                iso["MassUncertainty"] = unc
            elseif k == "Isotopic Composition"
                m = match(r"(.*?)\((.*?)#?\)", v)
                if m !== nothing
                    value = parse(Float64, m[1])
                    unc = parse(Float64, string(
                        replace(m[1], r"[0-9]" => "0")[1:end-length(m[2])], m[2]))
                elseif tryparse(Float64, v) !== nothing
                    value = parse(Float64, v)
                    unc = 0.0
                else
                    value = NaN
                    unc = NaN
                end
                iso["Composition"] = value
                iso["CompositionUncertainty"] = unc
            elseif k == "Standard Atomic Weight"
                m1 = match(r"\[(.*?),(.*?)\]", v)
                m2 = match(r"(.*?)\((.*?)\)", v)
                m3 = match(r"\[(.*?)\]", v)
                if m1 !== nothing
                    rcd["Weight"] = CONVENTIONAL_WEIGHT[rcd["Number"]]
                    rcd["WeightType"] = "Interval"
                    rcd["WeightLower"] = parse(Float64, m1[1])
                    rcd["WeightHigher"] = parse(Float64, m1[2])
                elseif m2 !== nothing
                    rcd["Weight"] = parse(Float64, m2[1])
                    unc = parse(Float64, string(
                        replace(m2[1], r"[0-9]" => "0")[1:end-length(m2[2])], m2[2]))
                    rcd["WeightType"] = "Uncertainty"
                    rcd["WeightUncertainty"] = unc
                elseif m3 !== nothing
                    rcd["Weight"] = parse(Int, m3[1])
                    rcd["WeightType"] = "MostStable"
                else
                    rcd["Weight"] = NaN
                    rcd["WeightType"] = "NotAvailable"
                end
            elseif k == "Notes"
                rcd["Notes"] = []
                v != " " && append!(rcd["Notes"], split(v, ","))
            end
            next = iterate(lines, state)
            next === nothing && break
        end
        revdict[rcd["Symbol"]] = rcd["Number"]
        if length(data) + 1 == rcd["Number"]
            if !isempty(data)
                e = data[end]
                # Monoisotopic mass
                comp = [
                    r["Composition"] for r in e["Isotopes"]
                    if r["Composition"] !== NaN]
                if !isempty(comp)
                    val, i = findmax(comp)
                    e["Monoisotopic"] = e["Isotopes"][i]["Mass"]
                    e["MonoisotopicUncertainty"
                        ] = e["Isotopes"][i]["MassUncertainty"]
                else
                    e["Monoisotopic"] = NaN
                    e["MonoisotopicUncertainty"] = NaN
                end
            end
            rcd["Isotopes"] = [iso]
            push!(data, rcd)
            # New entry
            if next === nothing
                e = data[end]
                # Monoisotopic mass
                comp = [
                    r["Composition"] for r in e["Isotopes"]
                    if r["Composition"] !== NaN]
                if !isempty(comp)
                    val, i = findmax(comp)
                    e["Monoisotopic"] = e["Isotopes"][i]["Mass"]
                    e["MonoisotopicUncertainty"
                        ] = e["Isotopes"][i]["MassUncertainty"]
                else
                    e["Monoisotopic"] = NaN
                    e["MonoisotopicUncertainty"] = NaN
                end
                break
            end
        else
            if rcd["Symbol"] != data[end]["Symbol"]
                iso["Symbol"] = rcd["Symbol"]
            end
            push!(data[end]["Isotopes"], iso)
        end
        next = iterate(lines, state)
    end
    println(length(data))
    YAML.write_file(dest, data)
    YAML.write_file(revdest, revdict)
end

run()
