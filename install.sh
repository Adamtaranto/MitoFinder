#!/bin/sh

wd=$(pwd)
chmod 764 bin/*/*
chmod 764 bin/*/*/*

# Might need to fix refs to this
#ln -s "$p"/spades.py "$p"/metaspades.py

cd "$wd"/mitfi
tar -xvf infernal-1.0.2.tar.gz
cd infernal-1.0.2
./configure  --prefix=$(pwd)/exec
make
make install

cd "$wd"/
touch install.sh.ok
