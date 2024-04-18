#!/bin/bash
cmake --build build
rm 2.txt
cd build
./txtreader BEM_INPUT_1_132273.txt
mv 2.txt ./..
