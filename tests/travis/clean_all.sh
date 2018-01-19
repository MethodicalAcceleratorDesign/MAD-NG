#! /bin/bash

# mad
make -C src -f Makefile.linux cleanall

# libs
make -C lib/luajit clean
make -C lib/lpeg   clean
make -C lib/lfs1   clean
make -C lib/lapack clean
make -C lib/fftw3  clean
make -C lib/nfft3  clean
make -C lib/nlopt2 clean

rm lib/{fftw3,nfft3,nlopt2}/Makefile

# rm -rf lib/{fftw3,nfft3,nlopt2,lapack,luajit}

# luarocks, etc
rm -rf install
