
using PackageCompiler

target_dir = "/opt/julia/moleculargraph"
# target_dir = replace(target_dir, "\\"=>"/")  # Change Windows paths to use "/"

println("Creating library in $target_dir")

PackageCompiler.create_library(
    ".", target_dir;
    lib_name="libmoleculargraph",
    precompile_execution_file=["./generate_precompile.jl"],
    header_files=["./libmoleculargraph.h"],
    incremental=true, force=true
)
