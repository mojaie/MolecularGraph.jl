JULIA ?= $(shell which julia)
JULIAC ?= $(shell $(JULIA) -e 'print(normpath(joinpath(Sys.BINDIR, Base.DATAROOTDIR, "julia", "juliac.jl")))')
OUTPUT := /opt/moleculargraphjl/libmoleculargraph.dylib
.PHONY: docs, clean, build

docs:
	@julia --project=docs/ --color=yes docs/make.jl

clean:
	rm -rf $(OUTPUT)

build: clean
	$(JULIA) --project=. $(JULIAC) --output-lib $(OUTPUT) --compile-ccallable src/MolecularGraph.jl

# Does not work for now
build-trim: clean
	$(JULIA) --project=. $(JULIAC) --experimental --output-lib $(OUTPUT) --trim=safe --compile-ccallable src/MolecularGraph.jl
