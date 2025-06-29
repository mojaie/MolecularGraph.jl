# Makefile

JULIA ?= julia
DLEXT := $(shell $(JULIA) --startup-file=no -e 'using Libdl; print(Libdl.dlext)')

DEST="/opt/julia/moleculargraph"

MYLIB_INCLUDES = $(DEST)/include/julia_init.h $(DEST)/include/libmoleculargraph.h
MYLIB_PATH := $(DEST)/lib/libmoleculargraph.$(DLEXT)

$(MYLIB_PATH) $(MYLIB_INCLUDES): ./build.jl ./libmoleculargraph.jl
	$(JULIA) --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate(); include("./libmoleculargraph.jl"); include("./build.jl")'


# For MacOS testing

link:
	install_name_tool -add_rpath $(DEST)/lib $(DEST)/lib/libmoleculargraph.$(DLEXT)
	install_name_tool -add_rpath $(DEST)/lib/julia $(DEST)/lib/libmoleculargraph.$(DLEXT)

unlink:
	install_name_tool -delete_rpath $(DEST)/lib $(DEST)/lib/libmoleculargraph.$(DLEXT)
	install_name_tool -delete_rpath $(DEST)/lib/julia $(DEST)/lib/libmoleculargraph.$(DLEXT)


.PHONY: clean
clean:
	$(RM) *~ *.o *.$(DLEXT)
	$(RM) -Rf $(TARGET)
