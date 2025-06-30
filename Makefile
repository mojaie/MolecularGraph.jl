JULIA ?= julia
JULIAC ?= $(shell $(JULIA) -e 'print(normpath(joinpath(Sys.BINDIR, Base.DATAROOTDIR, "julia", "juliac.jl")))')
DLEXT := $(shell $(JULIA) --startup-file=no -e 'using Libdl; print(Libdl.dlext)')
OUTPUT := /opt/julia/juliac-libmolgraphjl/libmolgraphjl.$(DLEXT)
.PHONY: docs, clean, build

docs:
	@julia --project=docs/ --color=yes docs/make.jl

clean:
	rm -rf $(OUTPUT)

build: clean
	$(JULIA) --project=./build $(JULIAC) --output-lib $(OUTPUT) --compile-ccallable ./build/src/LibMolGraphJL.jl

# Does not work for now
build-trim: clean
	$(JULIA) --project=. $(JULIAC) --experimental --output-lib $(OUTPUT) --trim=safe --compile-ccallable ./build/src/LibMolGraphJL.jl
