import PackageCompiler
using PackageCompiler: get_extra_linker_flags, julia_libdir, julia_private_libdir, ldlibs, bitflag, march, run_compiler
target_dir = "$(@__DIR__)/../compiled"
target_dir = replace(target_dir, "\\"=>"/")  # Change Windows paths to use "/"

# MacOS workaround
# https://github.com/JuliaLang/PackageCompiler.jl/issues/738

function PackageCompiler.create_sysimg_from_object_file(object_files::Vector{String},
        sysimage_path::String;
        version,
        compat_level::String,
        soname::Union{Nothing, String})

    if soname === nothing && (Sys.isunix() && !Sys.isapple())
        soname = basename(sysimage_path)
    end
    mkpath(dirname(sysimage_path))
    # Prevent compiler from stripping all symbols from the shared lib.
    o_file_flags = Sys.isapple() ? `-Wl,-all_load $object_files -Wl,-ld_classic` : `-Wl,--whole-archive $object_files -Wl,--no-whole-archive`
    extra = get_extra_linker_flags(version, compat_level, soname)
    cmd = `$(bitflag()) $(march()) -shared -L$(julia_libdir()) -L$(julia_private_libdir()) -o $sysimage_path $o_file_flags $(Base.shell_split(ldlibs())) $extra`
    run_compiler(cmd; cplusplus=true)
    return nothing
end

println("Creating library in $target_dir")
PackageCompiler.create_library(
    "../", target_dir;
    lib_name="libmoleculargraph",
    precompile_execution_file=["$(@__DIR__)/generate_precompile.jl"],
    header_files=["$(@__DIR__)/libmoleculargraph.h"],
    incremental=true
)
