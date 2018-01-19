#! /bin/bash
# NOTE: this must be built after luajit!
set -ex
HERE=$(readlink -f "$(dirname $BASH_SOURCE)")
name=lpeg-1.0.1
cd lib
[[ -e $name.tar.gz  ]] || wget http://www.inf.puc-rio.br/~roberto/lpeg/$name.tar.gz
[[ -e lpeg          ]] || (tar xzf $name.tar.gz && mv $name lpeg)
cd lpeg
patch -N -p0 <${HERE}/patch_lpeg.diff || true
make linux
cp liblpeg.a ../../bin/linux
cp re.lua ../../src/madl_regex.lua
