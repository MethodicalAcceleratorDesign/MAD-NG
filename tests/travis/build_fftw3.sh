#! /bin/bash
set -ex
cd lib
wget ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4.tar.gz
tar xzf fftw-3.3.4.tar.gz
mv fftw-3.3.4 fftw3
cd fftw3
./configure --disable-shared
make
cp .libs/libfftw3.a ../../bin/linux
