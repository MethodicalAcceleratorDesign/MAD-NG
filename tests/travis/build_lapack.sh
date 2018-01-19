#! /bin/bash
set -ex
HERE=$(readlink -f "$(dirname $BASH_SOURCE)")
[[ -e lib/lapack ]] || svn co https://icl.utk.edu/svn/lapack-dev/lapack/trunk lib/lapack
cd lib/lapack
cp make.inc.example make.inc
patch -p0 <${HERE}/patch_lapack.diff
make lapack_install lapacklib blaslib
cp liblapack.a librefblas.a ../../bin/linux
