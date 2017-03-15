#! /bin/bash
set -ex
git clone http://github.com/MethodicalAcceleratorDesign/LuaJIT.git lib/luajit
cd lib/luajit
git checkout mad-patch
make clean
make amalg PREFIX=`pwd`
make install PREFIX=`pwd`
mv bin/luajit{-2.1.0-beta2,}
cp src/libluajit.a ../../bin/linux
