# Makefile

JULIA ?= julia
DLEXT := $(shell $(JULIA) --startup-file=no -e 'using Libdl; print(Libdl.dlext)')

DEST="/opt/julia/libmolgraphjl"

MYLIB_INCLUDES = $(DEST)/include/julia_init.h $(DEST)/include/libmolgraphjl.h
MYLIB_PATH := $(DEST)/lib/libmolgraphjl.$(DLEXT)

$(MYLIB_PATH) $(MYLIB_INCLUDES): ./build.jl
	$(JULIA) --startup-file=no --project=. -e 'using Pkg; Pkg.develop(path=".."); Pkg.instantiate(); include("./build.jl")'
