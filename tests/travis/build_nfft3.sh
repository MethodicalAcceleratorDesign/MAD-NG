#! /bin/bash
set -ex
cd lib
wget http://www.nfft.org/download/nfft-3.3.1.tar.gz
tar xzf nfft-3.3.1.tar.gz
mv nfft-3.3.1 nfft3
cd nfft3
./configure --disable-shared \
            --with-fftw3=`pwd`/../fftw3 \
            --with-fftw3-libdir=`pwd`/../fftw3/.libs \
            --with-fftw3-includedir=`pwd`/../fftw3/api
make
cp .libs/libnfft3.a ../../bin/linux
