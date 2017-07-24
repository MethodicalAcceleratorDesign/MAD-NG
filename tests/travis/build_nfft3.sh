#! /bin/bash
set -ex
cd lib
[[ -e nfft-3.3.1.tar.gz ]] || wget http://www.nfft.org/download/nfft-3.3.1.tar.gz
[[ -e nfft3             ]] || (tar xzf nfft-3.3.1.tar.gz && \
                               mv nfft-3.3.1 nfft3)
cd nfft3
[[ -e Makefile          ]] || \
./configure --disable-shared --disable-examples --disable-applications \
            --with-fftw3=`pwd`/../fftw3 \
            --with-fftw3-libdir=`pwd`/../fftw3/.libs \
            --with-fftw3-includedir=`pwd`/../fftw3/api
make
cp .libs/libnfft3.a ../../bin/linux
