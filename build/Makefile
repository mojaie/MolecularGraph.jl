# Makefile

JULIA ?= julia
DLEXT := $(shell $(JULIA) --startup-file=no -e 'using Libdl; print(Libdl.dlext)')

TARGET="../compiled"
DEST="/opt/moleculargraphjl"

MYLIB_INCLUDES = $(TARGET)/include/julia_init.h $(TARGET)/include/libmoleculargraph.h
MYLIB_PATH := $(TARGET)/lib/libmoleculargraph.$(DLEXT)

$(MYLIB_PATH) $(MYLIB_INCLUDES): ./build.jl ../src/MolecularGraph.jl
	$(JULIA) --startup-file=no --project=../ -e 'using Pkg; Pkg.instantiate()'
	$(JULIA) --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate(); include("./build.jl")'


# For MacOS testing

link:
	mkdir -p $(DEST)
	cp -a $(TARGET)/* $(DEST)
	install_name_tool -add_rpath $(DEST)/lib $(DEST)/lib/libmoleculargraph.$(DLEXT)
	install_name_tool -add_rpath $(DEST)/lib/julia $(DEST)/lib/libmoleculargraph.$(DLEXT)

unlink:
	install_name_tool -delete_rpath $(DEST)/lib $(DEST)/lib/libmoleculargraph.$(DLEXT)
	install_name_tool -delete_rpath $(DEST)/lib/julia $(DEST)/lib/libmoleculargraph.$(DLEXT)


.PHONY: clean
clean:
	$(RM) *~ *.o *.$(DLEXT)
	$(RM) -Rf $(TARGET)
