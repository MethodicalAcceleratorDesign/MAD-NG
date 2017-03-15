#! /bin/bash
set -ex
svn co https://icl.utk.edu/svn/lapack-dev/lapack/trunk lib/lapack
cd lib/lapack
cp make.inc.example make.inc
make clean
make lapack_install lapacklib blaslib
cp liblapack.a librefblas.a ../../bin/linux
