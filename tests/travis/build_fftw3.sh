#! /bin/bash
set -ex
cd lib
name=fftw-3.3.6-pl2
[[ -e $name.tar.gz  ]] || wget ftp://ftp.fftw.org/pub/fftw/$name.tar.gz
[[ -e fftw3         ]] || (tar xzf $name.tar.gz && mv $name fftw3)
cd fftw3
[[ -e Makefile      ]] || ./configure --disable-shared
make
cp .libs/libfftw3.a ../../bin/linux
