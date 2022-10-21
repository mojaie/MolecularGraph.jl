
libmoleculargraph
===================================================

A library builder for calling some functions of MolecularGraph.jl from external applications. See [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl) for details.

### tested environment

```
julia> versioninfo()
Julia Version 1.8.2
Commit 36034abf260 (2022-09-29 15:21 UTC)
Platform Info:
  OS: macOS (x86_64-apple-darwin21.4.0)
  CPU: 4 × Intel(R) Core(TM) i5-8210Y CPU @ 1.60GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-13.0.1 (ORCJIT, skylake)
  Threads: 1 on 4 virtual cores
```

- ./scripts/test_lib.py worked with Python 3.10.0
- did not work with Julia version 1.7.3 due to a link error.