#!/bin/bash

rm *.o *.so *_wrap.c librainbow.py
# configuration file to wrap c rainbow library to python
# create wrap with swig
swig -python -py3 librainbow.i

# compile library
gcc -fpic -I/local/software/anaconda3/include/python3.5m -I/local/software/anaconda3/lib/python3.5/site-packages/numpy/core/include -I./include -c ./src/rainbow_read_raw.c ./src/rainbow_compress_raw.c ./src/qUncompress.c ./src/qCompress.c librainbow_wrap.c
gcc -shared -o _librainbow.so qUncompress.o qCompress.o rainbow_read_raw.o rainbow_compress_raw.o librainbow_wrap.o -lz