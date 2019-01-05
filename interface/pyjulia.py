
from julia import Main

code = '''
import Pkg
Pkg.activate(".")
using MolecularGraph
'''

Main.eval(code)


def sdfilereader(path):
    return Main.sdfilereader(path)


def sdftomol(path):
    return Main.sdftomol(path)


def drawsvg(mol, width, height):
    return Main.drawsvg_b(mol, width, height)
