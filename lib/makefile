#
# All Targets
#

all-lnx: lpeg-lnx
all-osx: lpeg-osx
all-win: lpeg-win

#
# LuaJIT (see README.luajit)
#

.PHONY: luajit

luajit-lnx: luajit
	rm -rf "$(PWD)/luajit/Linux" ; \
	make -C luajit install PREFIX="$(PWD)/luajit/Linux"

luajit-osx: luajit
	rm -rf "$(PWD)/luajit/OSX" ; \
	make -C luajit install PREFIX="$(PWD)/luajit/OSX"
 
luajit-win: luajit
	rm -rf "$(PWD)/luajit/Windows" ; \
	make -C luajit install PREFIX="$(PWD)/luajit/Windows"

luajit:
	make -C luajit clean ; \
	make -C luajit amalg

#
# LPEG (see README.lpeg)
#

lpeg-osx:
	make -C lpeg macosx LUADIR=../luajit/OSX/include/luajit-2.1 LIBNAME=lpeg-OSX.so

lpeg-lnx:
	make -C lpeg linux LUADIR=../luajit/Linux/include/luajit-2.1 LIBNAME=lpeg-Linux.so

lpeg-win:
	make -C lpeg windows LUADIR=../luajit/Linux/include/luajit-2.1 LIBNAME=lpeg-Windows.so
