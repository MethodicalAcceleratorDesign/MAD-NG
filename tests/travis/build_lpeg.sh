#! /bin/bash
# NOTE: this must be built after luajit!
set -ex
cd lib/lpeg
make clean
make
cp liblpeg.a ../../bin/linux
