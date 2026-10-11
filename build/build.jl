
using PackageCompiler
import PackageCompiler

target_dir = "/opt/julia/libmolgraphjl"
# target_dir = replace(target_dir, "\\"=>"/")  # Change Windows paths to use "/"

println("Creating library in $target_dir")

# TODO: PackageCompiler workaround for link error

function PackageCompiler.rpath_sysimage()
    Sys.iswindows() ? `` :
    Sys.isapple()   ? `-Wl,-rpath,'@loader_path' -Wl,-rpath,'@loader_path/julia'` :
                      `-Wl,-rpath,\$ORIGIN:\$ORIGIN/julia`
end


PackageCompiler.create_library(
    ".", target_dir;
    lib_name="libmolgraphjl",
    precompile_execution_file=["./generate_precompile.jl"],
    header_files=["./libmolgraphjl.h"],
    incremental=true, force=true
)
