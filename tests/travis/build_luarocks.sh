#! /bin/bash
set -ex

MAD_ROOT=$(readlink -f "$(dirname $BASH_SOURCE)"/../..)
LUA_ROOT=$MAD_ROOT/lib/luajit

cd lib

# acquire source
LUAROCKS=luarocks-2.4.2
[[ -e $LUAROCKS.tar.gz ]] || wget http://luarocks.org/releases/$LUAROCKS.tar.gz
[[ -e $LUAROCKS        ]] || tar -xzf $LUAROCKS.tar.gz

# build
cd $LUAROCKS
./configure --prefix=${PREFIX-$MAD_ROOT/install} \
            --with-lua-bin=$LUA_ROOT/bin \
            --with-lua-include=$LUA_ROOT/include/luajit-2.1 \
            --with-lua=$LUA_ROOT \
            --lua-suffix=jit
make clean
make build
make install
