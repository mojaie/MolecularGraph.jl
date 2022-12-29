using PackageCompiler

target_dir = "$(@__DIR__)/../compiled"
target_dir = replace(target_dir, "\\"=>"/")  # Change Windows paths to use "/"

println("Creating library in $target_dir")
PackageCompiler.create_library(
    "../", target_dir;
    lib_name="libmoleculargraph",
    precompile_execution_file=["$(@__DIR__)/generate_precompile.jl"],
    header_files=["$(@__DIR__)/libmoleculargraph.h"],
)
