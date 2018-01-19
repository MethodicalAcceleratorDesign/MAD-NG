#! /bin/bash

# Do this manually:
# cd src; make cleanall -f Makefile.linux
# git clean -xdf
# git stash

set -ex

./tests/travis/build_luajit.sh
./tests/travis/build_lpeg.sh
./tests/travis/build_lfs.sh

./tests/travis/build_lapack.sh
./tests/travis/build_fftw3.sh
./tests/travis/build_nfft3.sh
./tests/travis/build_nlopt2.sh

cd src
make -f Makefile.linux
cd ..

./tests/travis/build_luarocks.sh
. tests/travis/activate

ln -sf $(readlink -f src/mad) install/bin
ln -sf $(readlink -f lib/luajit/bin/luajit) install/bin/luajit
ln -sf luajit install/bin/lua
ln -sf luajit install/bin/lua51

luarocks install luacov
luarocks install cluacov
luarocks install lyaml
