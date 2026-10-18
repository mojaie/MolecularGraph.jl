# Python >= 3.10
#
# Test of libmolgraphjl built with JuliaC.
#
#   python3 python_test.py <bundle dir> [results.json]
#
# <bundle dir> is the directory given to `juliac --bundle` (contains lib/libmolgraphjl.*).
# If results.json is given, the cases and results are written to the file, which can be
# compared with the results of the Julia module by `reference.jl`.

import base64
from ctypes import CDLL, RTLD_GLOBAL, c_char_p, c_double, c_int, c_uint, c_void_p, string_at
import json
import math
from pathlib import Path
import platform
import sys

ASSETS_DIR = Path(__file__).resolve().parents[3] / "assets" / "test"

# Return types of the exported functions ("str": char* to be released by free_string)
FUNCTIONS = {
    "smilestomol": ([c_char_p, c_char_p], "str"),
    "smartstomol": ([c_char_p], "str"),
    "sdftomol": ([c_char_p, c_char_p], "str"),
    "vertex_count": ([c_char_p], "int"),
    "edge_count": ([c_char_p], "int"),
    "inchikey": ([c_char_p], "str"),
    "standard_weight": ([c_char_p], "double"),
    "molblock": ([c_char_p], "str"),
    "sdfmolblock": ([c_char_p], "str"),
    "drawsvg": ([c_char_p, c_char_p], "str"),
    "drawpng": ([c_char_p, c_uint, c_uint, c_char_p], "str"),
    "has_exact_match": ([c_char_p, c_char_p, c_char_p], "int"),
    "has_substruct_match": ([c_char_p, c_char_p, c_char_p], "int"),
    "tdmcis_size": ([c_char_p, c_char_p, c_char_p], "int"),
    "tdmces_size": ([c_char_p, c_char_p, c_char_p], "int"),
    "tdmcis_gls": ([c_char_p, c_char_p, c_char_p], "double"),
    "tdmces_gls": ([c_char_p, c_char_p, c_char_p], "double"),
    "last_error": ([], "str"),
}


class LibMolGraphJL:
    def __init__(self, bundle_dir):
        dlext = {"Darwin": "dylib", "Windows": "dll"}.get(platform.system(), "so")
        # The Julia runtime is initialized when the library is loaded (no jl_init needed)
        self.lib = CDLL(str(Path(bundle_dir) / "lib" / f"libmolgraphjl.{dlext}"), RTLD_GLOBAL)
        self.lib.free_string.argtypes = [c_void_p]
        self.lib.free_string.restype = None
        for name, (argtypes, ret) in FUNCTIONS.items():
            f = getattr(self.lib, name)
            f.argtypes = argtypes
            f.restype = {"str": c_void_p, "int": c_int, "double": c_double}[ret]

    def call(self, name, *args):
        args = [a.encode() if isinstance(a, str) else a for a in args]
        res = getattr(self.lib, name)(*args)
        if FUNCTIONS[name][1] == "str":
            s = string_at(res).decode("utf-8")
            self.lib.free_string(res)
            return s
        return res


def cases():
    # (name, function, args); {"ref": name} is replaced by the result of the case
    with open(ASSETS_DIR / "demo.mol") as f:
        demo = f.read()
    with open(ASSETS_DIR / "aspirin_3d.sdf") as f:
        aspirin3d = f.read()
    ref = lambda name: {"ref": name}
    return [
        ("aspirin", "smilestomol", ["CC(=O)OC1=CC=CC=C1C(=O)O", ""]),
        ("salicylic", "smilestomol", ["OC1=CC=CC=C1C(=O)O", "{}"]),
        ("salt", "smilestomol", ["CC(=O)[O-].[Na+]", '{"extract_largest_component": true}']),
        ("stereo1", "smilestomol", ["C[C@H](N)C(=O)O", ""]),
        ("stereo2", "smilestomol", ["C[C@@H](N)C(=O)O", ""]),
        ("invalid_smiles", "smilestomol", ["CC(C", ""]),
        ("invalid_smiles_error", "last_error", []),
        ("demo", "sdftomol", [demo, ""]),
        ("aspirin3d", "sdftomol", [aspirin3d, '{"remove_all_hydrogens": true}']),
        ("carboxyl", "smartstomol", ["[CX3](=O)[OX2H1]"]),
        ("rec", "smartstomol", ["[$(CO)]C"]),
        ("nested", "smartstomol", ["[$([CH]=[$(NOC)])]C=O"]),
        ("nestedor", "smartstomol", ["[$([CH]=[$(NOC),$(NO[Si])])]"]),
        ("nv", "vertex_count", [ref("aspirin")]),
        ("ne", "edge_count", [ref("aspirin")]),
        ("nv_salt", "vertex_count", [ref("salt")]),
        ("nv_invalid", "vertex_count", ["not a json"]),
        ("nv_invalid_error", "last_error", []),
        ("inchikey", "inchikey", [ref("aspirin")]),
        ("inchikey_demo", "inchikey", [ref("demo")]),
        ("weight", "standard_weight", [ref("aspirin")]),
        ("weight_invalid", "standard_weight", ["{}"]),
        ("molblock", "molblock", [ref("aspirin")]),
        ("sdfmolblock", "sdfmolblock", [ref("demo")]),
        ("svg", "drawsvg", [ref("aspirin"), ""]),
        ("svg_demo", "drawsvg", [ref("demo"), '{"show_carbon": "all", "bgcolor": "#F0F0F0", "viewbox": false}']),
        ("svg_query", "drawsvg", [ref("carboxyl"), ""]),  # not supported
        ("svg_query_error", "last_error", []),
        ("png", "drawpng", [ref("aspirin"), 200, 200, ""]),
        ("png_demo", "drawpng", [ref("demo"), 300, 200, '{"bgopacity": 0.5}']),
        ("exact_self", "has_exact_match", [ref("aspirin"), ref("aspirin"), ""]),
        ("exact_stereo", "has_exact_match", [ref("stereo1"), ref("stereo2"), '{"stereo": true}']),
        ("exact_nostereo", "has_exact_match", [ref("stereo1"), ref("stereo2"), "{}"]),
        ("substr", "has_substruct_match", [ref("aspirin"), ref("salicylic"), ""]),
        ("substr_rev", "has_substruct_match", [ref("salicylic"), ref("aspirin"), ""]),
        ("substr_query", "has_substruct_match", [ref("aspirin"), ref("carboxyl"), ""]),
        ("substr_rec", "has_substruct_match", [ref("aspirin"), ref("rec"), '{"timeout": 10}']),
        ("query_query", "has_substruct_match", [ref("nested"), ref("nestedor"), ""]),
        ("query_query_rev", "has_substruct_match", [ref("nestedor"), ref("nested"), ""]),
        ("query_first", "has_substruct_match", [ref("carboxyl"), ref("aspirin"), ""]),
        ("query_first_error", "last_error", []),
        ("mcis", "tdmcis_size", [ref("aspirin"), ref("salicylic"), ""]),
        ("mces", "tdmces_size", [ref("aspirin"), ref("salicylic"), ""]),
        ("mces_diam", "tdmces_size", [ref("aspirin"), ref("demo"), '{"diameter": 4, "tolerance": 1}']),
        ("mcis_gls", "tdmcis_gls", [ref("aspirin"), ref("salicylic"), ""]),
        ("mces_gls", "tdmces_gls", [ref("aspirin"), ref("salicylic"), ""]),
        ("mcs_query", "tdmces_size", [ref("aspirin"), ref("carboxyl"), ""]),
    ]


def run(lib):
    results = {}
    resolved_cases = []
    for name, func, args in cases():
        args = [results[a["ref"]] if isinstance(a, dict) else a for a in args]
        results[name] = lib.call(func, *args)
        resolved_cases.append({"name": name, "func": func, "args": args})
    return resolved_cases, results


def check(results):
    # Sanity checks independent of the Julia results
    r = results
    assert json.loads(r["aspirin"])["graph"], "aspirin should have edges"
    assert "error_smiles" in json.loads(r["invalid_smiles"])["gprops"]["logs"]
    assert r["invalid_smiles_error"] != ""
    assert (r["nv"], r["ne"], r["nv_salt"]) == (13, 13, 4)
    assert r["nv_invalid"] == -1 and r["nv_invalid_error"] != ""
    assert r["inchikey"] == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
    assert abs(r["weight"] - 180.16) < 1e-6
    assert r["molblock"].rstrip().endswith("M  END")
    assert r["sdfmolblock"].rstrip().endswith("$$$$")
    assert r["svg"].startswith("<svg")
    assert r["svg_query"] == "" and r["svg_query_error"] != ""
    png = base64.b64decode(r["png"])
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    assert (r["exact_self"], r["exact_stereo"], r["exact_nostereo"]) == (1, 0, 1)
    assert (r["substr"], r["substr_rev"], r["substr_query"], r["substr_rec"]) == (1, 0, 1, 1)
    assert (r["query_query"], r["query_query_rev"]) == (1, 0)
    assert r["query_first"] == -1 and r["query_first_error"] != ""
    assert r["mcs_query"] == -1
    assert 0 < r["mcis_gls"] <= 1 and 0 < r["mces_gls"] <= 1


if __name__ == "__main__":
    lib = LibMolGraphJL(sys.argv[1])
    resolved_cases, results = run(lib)
    for name, res in results.items():
        s = repr(res)
        print(f"{name}: {s[:80]}{'...' if len(s) > 80 else ''}")
    check(results)
    print("sanity checks passed")
    if len(sys.argv) > 2:
        dump = {k: (None if isinstance(v, float) and math.isnan(v) else v) for k, v in results.items()}
        with open(sys.argv[2], "w") as f:
            json.dump({"cases": resolved_cases, "results": dump}, f)
