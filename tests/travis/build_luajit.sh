#! /bin/bash
set -ex
[[ -e lib/luajit ]] || git clone http://github.com/MethodicalAcceleratorDesign/LuaJIT.git lib/luajit
cd lib/luajit
git checkout mad-patch
git pull
make amalg PREFIX=`pwd`
make install PREFIX=`pwd`
mv bin/luajit{-2.1.0-beta3,}
cp src/libluajit.a ../../bin/linux
