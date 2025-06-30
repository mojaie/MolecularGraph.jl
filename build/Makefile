# Makefile

JULIA ?= julia
DLEXT := $(shell $(JULIA) --startup-file=no -e 'using Libdl; print(Libdl.dlext)')

DEST="/opt/julia/libmolgraphjl"

MYLIB_INCLUDES = $(DEST)/include/julia_init.h $(DEST)/include/libmolgraphjl.h
MYLIB_PATH := $(DEST)/lib/libmolgraphjl.$(DLEXT)

$(MYLIB_PATH) $(MYLIB_INCLUDES): ./build.jl
	$(JULIA) --startup-file=no --project=. -e 'using Pkg; Pkg.develop(path=".."); Pkg.instantiate(); include("./build.jl")'


# For MacOS testing

link:
	install_name_tool -add_rpath $(DEST)/lib $(DEST)/lib/libmolgraphjl.$(DLEXT)
	install_name_tool -add_rpath $(DEST)/lib/julia $(DEST)/lib/libmolgraphjl.$(DLEXT)

unlink:
	install_name_tool -delete_rpath $(DEST)/lib $(DEST)/lib/libmolgraphjl.$(DLEXT)
	install_name_tool -delete_rpath $(DEST)/lib/julia $(DEST)/lib/libmolgraphjl.$(DLEXT)


.PHONY: clean
clean:
	$(RM) *~ *.o *.$(DLEXT)
	$(RM) -Rf $(TARGET)
