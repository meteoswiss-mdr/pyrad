#!/bin/bash

# script to generate the MeteoSwiss Pyart library reference and remove intermediate outputs
# fvj 20.10.2016

# library reference for developers
rm -f ../pyart-mch_library_reference_dev.pdf
cd library_reference_developers

# make pdf
# make clean
# make latexpdf
# cp build/latex/*.pdf ../../.
# rm -rf build/*

# make html for github pages
# make clean
# make html
# rm -r ../../../docs/*
# mv build/html/* ../../../docs
# rm -rf build/*
cd ..

# library reference for users
rm -f ../pyart-mch_library_reference_users.pdf
cd library_reference_users

# # make pdf
# make clean
# make latexpdf
# cp build/latex/*.pdf ../../.
# rm -rf build/*
# cd ..

# make html for github pages
make clean
make html
git rm -r ../../../src/pyart/docs/*
mv build/html/* ../../../src/pyart/docs
rm -rf build/*
touch ../../../src/pyart/docs/.nojekyll
cd ..