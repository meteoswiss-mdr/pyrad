#!/bin/bash

# script to generate the pyrad library reference and remove intermediate outputs
# fvj 20.10.2016

# library reference for developers
rm -f ../pyrad_library_reference_dev.pdf
cd library_reference_developers
make clean

# make pdf
# make latexpdf
# cp build/latex/*.pdf ../../.
# rm -rf build/*
# cd ..
make html
rm -r ../../../docs/*
mv build/html/* ../../../docs

# # library reference for users
# rm -f ../pyrad_library_reference_users.pdf
# cd library_reference_users
# make clean
# make latexpdf
# cp build/latex/*.pdf ../../.
# rm -rf build/*
# cd ..