
import json
from pathlib import Path
from ctypes import CDLL, RTLD_GLOBAL, c_char_p, c_double, c_int


scriptdir = Path(__file__).parent

julia = CDLL(scriptdir.parent / "compiled/lib/libmoleculargraph.dylib", RTLD_GLOBAL)
try:
    julia.jl_init_with_image__threading(
        bytes(scriptdir.parent / "compiled/lib"),
        b"libmoleculargraph.dylib"
    )
except AttributeError:
    julia.jl_init_with_image(
        bytes(scriptdir.parent / "compiled/lib"),
        b"libmoleculargraph.dylib"
    )

julia.smilestomol.argtypes = [c_char_p]
julia.smilestomol.restype = c_char_p
julia.inchikey.argtypes = [c_char_p]
julia.inchikey.restype = c_char_p
smiles = b"CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
mol1 = julia.smilestomol(smiles)
ikey = julia.inchikey(mol1)
print(ikey.decode("utf-8"))

julia.sdftomol.argtypes = [c_char_p]
julia.sdftomol.restype = c_char_p
julia.standard_weight.argtypes = [c_char_p]
julia.standard_weight.restype = c_double

with open(scriptdir.parent / "assets/test/demo.mol") as f:
    sdf = f.read()
mol2 = julia.sdftomol(sdf.encode())
mw = julia.standard_weight(mol2)
print(mw)

julia.has_exact_match.argtypes = [c_char_p, c_char_p, c_char_p]
julia.has_exact_match.restype = c_int
julia.has_substruct_match.argtypes = [c_char_p, c_char_p, c_char_p]
julia.has_substruct_match.restype = c_int
julia.tdmcis.argtypes = [c_char_p, c_char_p, c_char_p]
julia.tdmcis.restype = c_int
julia.tdmces.argtypes = [c_char_p, c_char_p, c_char_p]
julia.tdmces.restype = c_int

match1 = julia.has_exact_match(mol1, mol2, json.dumps({}).encode())
print(match1)
match2 = julia.has_substruct_match(mol1, mol2, json.dumps({}).encode())
print(match2)
cnt1 = julia.tdmcis(mol1, mol2, json.dumps({"diameter": 2, "tolerance": 1}).encode())
print(cnt1)
cnt2 = julia.tdmces(mol1, mol2, json.dumps({"diameter": 2, "tolerance": 1}).encode())
print(cnt2)

julia.jl_atexit_hook(0)
