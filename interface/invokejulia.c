
#include <stdio.h>
#include <julia.h>
#include "invokejulia.h"

JULIA_DEFINE_FAST_TLS()


jl_value_t * smilestomol(char *smiles)
{
    jl_init();
    jl_eval_string("using GraphMol");
    char buffer[1000];
    sprintf(buffer, "smilestomol(\"%s\")", smiles);
    jl_value_t *res = jl_eval_string(buffer);
    jl_atexit_hook(0);
    return res;
}

/*
$JULIA/share/julia/julia-config.jl --allflags

gcc -o invokejulia.so invokejulia.c -shared -I'/Applications/Julia-1.0.app/Contents/Resources/julia/include/julia' -DJULIA_ENABLE_THREADING=1 -fPIC -L'/Applications/Julia-1.0.app/Contents/Resources/julia/lib' -Wl,-rpath,'/Applications/Julia-1.0.app/Contents/Resources/julia/lib' -Wl,-rpath,'/Applications/Julia-1.0.app/Contents/Resources/julia/lib/julia' -ljulia
*/
