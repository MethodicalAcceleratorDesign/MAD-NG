#! /bin/bash
set -ex
[[ -e lib/lfs1 ]] || git clone \
    https://github.com/MethodicalAcceleratorDesign/luafilesystem.git lib/lfs1
cd lib/lfs1
make lfs.a
cp liblfs.a ../../bin/linux
