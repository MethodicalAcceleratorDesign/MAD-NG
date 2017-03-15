#! /bin/bash
set -ex
cd lib
[[ -e fftw-3.3.4.tar.gz ]] || wget ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4.tar.gz
[[ -e fftw3             ]] || (tar xzf fftw-3.3.4.tar.gz && \
                               mv fftw-3.3.4 fftw3)
cd fftw3
./configure --disable-shared
make clean
make
cp .libs/libfftw3.a ../../bin/linux
