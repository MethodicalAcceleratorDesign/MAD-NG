#! /bin/bash
set -ex
cd lib
name=nfft-3.3.2
[[ -e $name.tar.gz  ]] || wget http://www.nfft.org/download/$name.tar.gz
[[ -e nfft3         ]] || (tar xzf $name.tar.gz && mv $name nfft3)
cd nfft3
[[ -e Makefile      ]] || \
./configure --disable-shared --disable-examples --disable-applications \
            --with-fftw3=`pwd`/../fftw3 \
            --with-fftw3-libdir=`pwd`/../fftw3/.libs \
            --with-fftw3-includedir=`pwd`/../fftw3/api
make
cp .libs/libnfft3.a ../../bin/linux
