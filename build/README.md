
libmoleculargraph
===================================================

A library builder for calling some functions of MolecularGraph.jl from external applications. See [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl) for details.


### Dockerfile (recommended)

- See `./docker-example`
- A test script `../scripts/python_test.py` worked with Python 3.10.0


### Other platforms

Building in environments shown below and other platforms are not supported.


#### For Intel Mac

- Building MacOS binary has some potential issues.
- Run `make` to build binary and then `sudo make link` to copy them to `/usr/local`
- Testing environment:

  ```
  julia> versioninfo()
  Julia Version 1.10.0
  Commit 3120989f39b (2023-12-25 18:01 UTC)
  Build Info:
    Official https://julialang.org/ release
  Platform Info:
    OS: macOS (x86_64-apple-darwin22.4.0)
    CPU: 4 × Intel(R) Core(TM) i5-8210Y CPU @ 1.60GHz
    WORD_SIZE: 64
    LIBM: libopenlibm
    LLVM: libLLVM-15.0.7 (ORCJIT, skylake)
    Threads: 1 on 4 virtual cores
  ```

  - did not work with Julia version 1.7.3 (and below) due to a link error.