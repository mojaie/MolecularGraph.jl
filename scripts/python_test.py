from ctypes import CDLL, RTLD_GLOBAL, POINTER, cast, c_char_p, c_ubyte, c_double, c_int, c_uint
import json
from pathlib import Path
import platform


def jl_init():
    dlext = "dylib" if platform.system() == "Darwin" else "so"
    libdir = Path("/usr/local/lib")
    jl = CDLL(libdir / f"libmoleculargraph.{dlext}", RTLD_GLOBAL)
    jl.jl_init_with_image(bytes(libdir), f"libmoleculargraph.{dlext}".encode())
    return jl


def jl_exit(jl):
    jl.jl_atexit_hook(0)


def smiles_to_mol(jl, smiles: str) -> bytes:
    jl.smilestomol.argtypes = [c_char_p]
    jl.smilestomol.restype = c_char_p
    return jl.smilestomol(smiles.encode(), r"{}".encode())


def smarts_to_mol(jl, smarts: str) -> bytes:
    jl.smartstomol.argtypes = [c_char_p]
    jl.smartstomol.restype = c_char_p
    return jl.smartstomol(smarts.encode())


def sdf_to_mol(jl, sdf: str) -> bytes:
    # sdf: SDFile string (open as f -> f.read())
    jl.sdftomol.argtypes = [c_char_p]
    jl.sdftomol.restype = c_char_p
    return jl.sdftomol(sdf.encode(), r"{}".encode())


def vertex_count(jl, mol: bytes) -> int:
    # mol: json.dumps(mol_dict).encode()
    jl.vertex_count.argtypes = [c_char_p]
    jl.vertex_count.restype = c_int
    return jl.vertex_count(mol)


def edge_count(jl, mol: bytes) -> int:
    # mol: json.dumps(mol_dict).encode()
    jl.edge_count.argtypes = [c_char_p]
    jl.edge_count.restype = c_int
    return jl.edge_count(mol)


def inchikey(jl, mol: bytes) -> str:
    # mol: json.dumps(mol_dict).encode()
    jl.inchikey.argtypes = [c_char_p]
    jl.inchikey.restype = c_char_p
    return jl.inchikey(mol).decode("utf-8")


def standard_weight(jl, mol: bytes) -> float:
    # mol: json.dumps(mol_dict).encode()
    jl.standard_weight.argtypes = [c_char_p]
    jl.standard_weight.restype = c_double
    return jl.standard_weight(mol)


def mol_to_svg(jl, mol: bytes) -> str:
    # mol: json.dumps(mol_dict).encode()
    jl.drawsvg.argtypes = [c_char_p]
    jl.drawsvg.restype = c_char_p
    return jl.drawsvg(mol).decode("utf-8")


def sdf_to_svg(jl, sdf: str) -> str:
    mol = sdf_to_mol(jl, sdf)
    return mol_to_svg(jl, mol)


def smiles_to_svg(jl, smiles: str) -> str:
    mol = smiles_to_mol(jl, smiles)
    return mol_to_svg(jl, mol)


def mol_to_png(jl, mol: bytes, width: int, height: int) -> bytes:
    # mol: json.dumps(mol_dict).encode()
    buf = (c_ubyte * 200_000).from_buffer(bytearray(200_000))  # 200 kB buffer, tentative
    jl.drawpng.argtypes = [c_ubyte * 200_000, c_char_p, c_uint, c_uint]
    jl.drawpng.restype = c_int
    size = jl.drawpng(buf, mol, c_uint(width), c_uint(height))
    return bytes(buf[:size])


def sdf_to_png(jl, sdf: str, width: int, height: int) -> bytes:
    mol = sdf_to_mol(jl, sdf)
    return mol_to_png(jl, mol, width, height)


def smiles_to_png(jl, smiles: str, width: int, height: int) -> bytes:
    mol = smiles_to_mol(jl, smiles)
    return mol_to_png(jl, mol, width, height)


def has_exact_match(jl, mol1: bytes, mol2: bytes, options: bytes) -> bool:
    # mol1, mol2, options: json.dumps(dict).encode()
    jl.has_exact_match.argtypes = [c_char_p, c_char_p, c_char_p]
    jl.has_exact_match.restype = c_int
    return jl.has_exact_match(mol1, mol2, options)


def has_substruct_match(jl, mol1: bytes, mol2: bytes, options: bytes) -> bool:
    # mol1, mol2, options: json.dumps(dict).encode()
    jl.has_substruct_match.argtypes = [c_char_p, c_char_p, c_char_p]
    jl.has_substruct_match.restype = c_int
    return jl.has_substruct_match(mol1, mol2, options)


def tdmcis_size(jl, mol1: bytes, mol2: bytes, options: bytes) -> int:
    # mol1, mol2, options: json.dumps(dict).encode()
    jl.tdmcis_size.argtypes = [c_char_p, c_char_p, c_char_p]
    jl.tdmcis_size.restype = c_int
    return jl.tdmcis_size(mol1, mol2, options)


def tdmces_size(jl, mol1: bytes, mol2: bytes, options: bytes) -> int:
    # mol1, mol2, options: json.dumps(dict).encode()
    jl.tdmces_size.argtypes = [c_char_p, c_char_p, c_char_p]
    jl.tdmces_size.restype = c_int
    return jl.tdmces_size(mol1, mol2, options)


if __name__ == "__main__":
    jl = jl_init()
    mol1 = smiles_to_mol(jl, "CC(=O)OC1=CC=CC=C1C(=O)O")
    mol2 = smiles_to_mol(jl, "OC1=CC=CC=C1C(=O)O")
    print(standard_weight(jl, mol1))
    print(inchikey(jl, mol1))
    print(has_substruct_match(jl, mol1, mol2, json.dumps({}).encode()))
    print(tdmces_size(jl, mol1, mol2, json.dumps({}).encode()))
    print(mol_to_svg(jl, mol1))
    jl_exit(jl)