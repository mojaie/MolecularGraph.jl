
Build package binary
===================================================

1. Test environment

    - MolecularGraph.jl 0.20.2
    - Julia 1.11.5
    - Python 3.13.2
    - MacOS 15.5 (M3, for local test)
    - Ubuntu 22.04.5 LTS (production)
      - Docker Engine 27.3.1

1. Create destination directory

    - e.g. `sudo mkdir /opt/julia` and `sudo chown <user>:<group> /opt/julia`

1. Build binary

    - PackageCompiler.jl

      See [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl) for details.

      ```sh
      cd $YOUR_WORKSPACE/MolecularGraph.jl
      cd build
      make
      ```

    - JuliaC.jl (experimental, see "Build with JuliaC" below)

1. Add library path

    add to your profile file (e.g. `.zprofile`)  
    - `export DYLD_FALLBACK_LIBRARY_PATH='/opt/julia/libmolgraphjl/lib:/opt/julia/libmolgraphjl/lib/julia'`

1. Test Python script

    ```sh
    cd $YOUR_WORKSPACE/MolecularGraph.jl
    python ./scripts/python_test.py
    ```

1. Create Docker image

    - See `./build/docker-example`


Build with JuliaC (experimental)
===================================================

`build/juliac` is a package for building a small shared library with
[JuliaC.jl](https://github.com/JuliaLang/JuliaC.jl) and `--trim=safe` (tested with Julia 1.13.2,
JuliaC 0.3.12). It provides the same functions as the PackageCompiler version (`build/src`)
with the following differences.

- The Julia runtime is initialized when the library is loaded. Do not call `jl_init_with_image`
  or `jl_atexit_hook`.
- `jl_init`/precompile statements (`generate_precompile.jl`) are not needed. All code reachable
  from the `@ccallable` functions is compiled ahead of time.
- Strings returned by the library are allocated by `malloc`. Release them with
  `free_string(char*)` (use `c_void_p` as `restype` in ctypes).
- On error, functions return a sentinel value (`-1` for integers, `NaN` for floats, an empty
  string, or an empty molecule with the message in `gprops.logs`), and the message can be
  obtained by `last_error()`.
- `drawpng` requires `Cairo` (included in the dependencies).

1. Install JuliaC

    ```sh
    julia -e 'using Pkg; Pkg.Apps.add("JuliaC")'
    ```

    `juliac` is installed in `~/.julia/bin`.

1. Build

    ```sh
    cd $YOUR_WORKSPACE/MolecularGraph.jl
    julia --project=build/juliac -e 'using Pkg; Pkg.instantiate()'
    juliac --output-lib libmolgraphjl --bundle /opt/julia/libmolgraphjl --trim=safe --experimental build/juliac
    ```

    The library is created at `/opt/julia/libmolgraphjl/lib/libmolgraphjl.{dylib,so}` with the
    Julia runtime and the binary dependencies (`lib/julia`). Unresolved calls are reported by
    the trim verifier and the build fails; see `.agents/decisions.md` for the patterns that
    need to be avoided.

1. Test

    ```sh
    python3 build/juliac/test/python_test.py /opt/julia/libmolgraphjl results.json
    julia --project=build/juliac build/juliac/test/reference.jl results.json
    ```

    `python_test.py` calls the library and runs sanity checks. `reference.jl` runs the same
    cases with the `LibMolGraphJL` module on Julia and checks that the results are identical.
