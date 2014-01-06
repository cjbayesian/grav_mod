#!/bin/bash

## Hack makefile ##

loc=`uname -n`
if [ $loc == "this-UX31E" ]
then
   g++ -w -I/home/cchivers/SchoolBackUp/c++/library/ -I/usr/share/R/include -o "$1" "src/$1.cpp" /home/cchivers/SchoolBackUp/c++/library/libCJ -lRmath -lm -fopenmp -O3 -ftree-vectorize -msse2
fi

if [ $loc == "chivers" ]
then
   g++ -w -I/home/cchivers/c++/library/ -I/usr/share/R/include -o "$1" "src/$1.cpp" /home/cchivers/c++/library/libCJ -lRmath -lm -fopenmp -O3 -ftree-vectorize -msse2
fi

exit 0
