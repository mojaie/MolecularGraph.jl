
from ctypes import CDLL, c_char_p

JL = CDLL("invokejulia.so")

JL.smilestomol.argtypes = [c_char_p]
print(JL.smilestomol("CCC(=O)CCO".encode()))

"""
def mol_to_svg(mol):
    JL.smilestomol.argtypes = None  # Pointer to MolGraph
    return JL.moltosvg(mol)
"""
