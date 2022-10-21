
from pathlib import Path
from ctypes import CDLL, RTLD_GLOBAL, c_char_p, c_double


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

julia.smilesmoldata.argtypes = [c_char_p]
julia.smilesmoldata.restype = c_char_p
teststr = "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"

res = julia.smilesmoldata(b"CCC")
print(res.decode("utf-8"))

julia.jl_atexit_hook(0)
