#
# Packages
#

packages := luajit lpeg

#
# Targets
#

.PHONY: $(packages)

all:    $(packages)

###############################################################################

#
# LuaJIT (see README.luajit)
#

luajit:
	make -C $@ amalg ; \
	make -C $@ install PREFIX="$(PWD)/$@/$(OS)"

#
# LPEG (see README.lpeg)
#

lpeg:
	make -C $@ $(TG) LUADIR=../luajit/$(OS)/include/luajit-2.1 LIBNAME=lpeg-$(OS).so

###############################################################################

#
# Cleaning
#

clean:
	$(foreach pkg,$(packages),make -C $(pkg) clean ;)

cleanall: clean
	rm -rf $(addprefix $(PWD)/luajit/,$(OSLIST))

#
# Detect platform (same as luajit jit.os)
#

ifeq ($(or $(windir),$(WINDIR)),)
OS := $(subst Darwin,OSX,$(shell uname -s))
else
OS := Windows
endif

OSLIST := Linux OSX Windows

ifeq ($(filter $(OS),$(OSLIST)),)
$(error Unsupported platform $(OS))
endif

TG = $(target-$(OS))

target-Linux   := linux
target-OSX     := macosx
target-Windows := windows
