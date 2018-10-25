#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Atom,
    getnumber,
    getname,
    getcolor,
    getweight,
    sethydrogen!,
    chargesign,
    markup,
    html


const PERIODIC_TABLE = YAML.load(open(
    joinpath(dirname(@__FILE__), "..", "..", "assets", "const", "periodictable.yaml")
))
const H_WEIGHT = PERIODIC_TABLE["H"]["std_weight"]


mutable struct Atom <: Node
    index::UInt16

    symbol::String
    charge::Int8
    multiplicity::UInt8
    mass::Union{UInt8, Nothing}

    Hcount::UInt8
    pi::UInt8
    aromatic::Bool
    Hdonor::Bool
    Hacceptor::Bool
    carbonylC::Bool
    lonepair::Bool

    wctype::UInt8
    patty::UInt8
    smiles_aromatic::Bool
    smiles_stereo::String
    coords::SVector{3, Float32}
    visible::Bool

    function Atom(symbol::AbstractString)
        atom = new()
        if symbol ∉ keys(PERIODIC_TABLE)
            throw(DescriptorError("Atom '$(symbol)' not supported"))
        end
        atom.symbol = symbol
        atom.charge = 0
        atom.multiplicity = 1
        atom.mass = nothing

        atom.Hcount = 0
        atom.pi = 0
        atom.aromatic = false
        atom.Hdonor = false
        atom.Hacceptor = symbol in ("N", "O", "F")
        atom.carbonylC = false
        atom.lonepair = false

        atom.smiles_aromatic = false
        atom.smiles_stereo = ""
        atom.visible = symbol != "C"
        atom
    end
end

Atom() = Atom("C")


function getnumber(atom::Atom)
    PERIODIC_TABLE[atom.symbol]["number"]
end


function getname(atom::Atom)
    PERIODIC_TABLE[atom.symbol]["name"]
end


function getcolor(atom::Atom)
    attr = PERIODIC_TABLE[atom.symbol]
    tuple(get(attr, "color", [0, 192, 192])...)
end


function getweight(atom::Atom)
    m = PERIODIC_TABLE[atom.symbol]["std_weight"]
    m + H_WEIGHT * atom.Hcount
end


function sethydrogen!(atom::Atom, Hs::UInt8)
    atom.Hcount = Hs
    atom.Hdonor = Hs > 0 && atom.symbol in ("N", "O")
end


function chargesign(atom::Atom)
    if atom.charge == 0
        return ""
    end
    sign = atom.charge > 0 ? "+" : "–" # en dash, not hyphen-minus
    num = abs(atom.charge)
    num > 1 ? "$(num)$(sign)" : sign
end


function markup(atom::Atom, direction::Symbol,
                substart::String, subend::String,
                supstart::String, supend::String)
    if atom.Hcount == 1
        text = "H"
    elseif atom.Hcount > 1
        text = string("H", substart, atom.Hcount, subend)
    else
        text = ""
    end
    chg = atom.charge == 0 ? "" : string(supstart, chargesign(atom), supend)
    seq = [atom.symbol, text, chg]
    if direction == :left
        seq = reverse(seq)
    end
    join(seq, "")
end


function html(atom::Atom, direction::Symbol)
    markup(atom, direction, "<sub>", "</sub>", "<sup>", "</sup>")
end
