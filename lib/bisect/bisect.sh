make clean
make amalg PREFIX=`pwd`
make install PREFIX=`pwd`
mv bin/luajit-2.1.0-beta2 bin/luajit
bin/luajit $1