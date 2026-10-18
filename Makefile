
OUTDIR := /opt/julia/libmolgraphjl
.PHONY: docs clean build

docs:
	julia --project=docs --color=yes docs/make.jl

clean:
	rm -rf $(OUTDIR) && mkdir $(OUTDIR)

build: clean
	julia --project=build/juliac -e 'using Pkg; Pkg.instantiate()' \
	&& juliac --output-lib libmolgraphjl --bundle $(OUTDIR) --trim=safe --experimental build/juliac

# Add Cairo and RDKitMinimalLib to the default env in advance
dtest:
	julia --project=. ext/tests_draw.jl
