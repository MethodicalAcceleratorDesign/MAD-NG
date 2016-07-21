make clean
make amalg PREFIX=`pwd`
make install PREFIX=`pwd`
cp bin/luajit-2.1.0-beta2 bin/luajit
bin/luajit $@