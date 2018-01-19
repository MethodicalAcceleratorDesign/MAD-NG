#! /bin/bash
set -ex
cd lib
name=nlopt-2.4.2
[[ -e $name.tar.gz  ]] || wget http://ab-initio.mit.edu/nlopt/$name.tar.gz
[[ -e nlopt2        ]] || (tar xzf $name.tar.gz && mv $name nlopt2)
cd nlopt2
[[ -e Makefile      ]] || ./configure --disable-shared
make
cp .libs/libnlopt.a ../../bin/linux/libnlopt2.a
