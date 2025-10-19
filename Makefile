JULIA ?= julia
JULIAC ?= $(shell $(JULIA) -e 'print(normpath(joinpath(Sys.BINDIR, Base.DATAROOTDIR, "julia", "juliac.jl")))')
DLEXT := $(shell $(JULIA) --startup-file=no -e 'using Libdl; print(Libdl.dlext)')
OUTDIR := /opt/julia/juliac-libmolgraphjl
OUTPUT := $(OUTDIR)/libmolgraphjl.$(DLEXT)
.PHONY: docs, clean, build

docs:
	@julia --project=docs/ --color=yes docs/make.jl

clean:
	rm -rf $(OUTDIR) && mkdir $(OUTDIR)

build: clean
	$(JULIA) --project=./build $(JULIAC) --output-lib $(OUTPUT) --compile-ccallable ./build/src/LibMolGraphJL.jl

# Does not work for now
build-trim: clean
	$(JULIA) --project=. $(JULIAC) --experimental --output-lib $(OUTPUT) --trim=safe --compile-ccallable ./build/src/LibMolGraphJL.jl

# Add Cairo and RDKitMinimalLib to the default env in advance
dtest:
	julia --project=. ext/tests_draw.jl

xtest:
	julia --project=. ext/tests_ext.jl