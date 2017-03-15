#! /bin/bash
# NOTE: this must be built after luajit!
set -ex
cd lib/lpeg
make
cp liblpeg.a ../../bin/linux
