#!/bin/sh

git clone git@github.com:pezmaster31/bamtools.git
mkdir build
cd build
cmake ..
make -j 6
sudo make install
sudo mv /usr/local/include/bamtools/* /usr/local/include/
sudo mv /usr/local/lib/bamtools/libbamtools.* /usr/local/lib
sudo ldconfig
