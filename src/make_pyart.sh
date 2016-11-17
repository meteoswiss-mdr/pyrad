#!/bin/bash

# script to build the MeteoSwiss Pyart
# fvj 17.11.2016

# remove previous built
echo 'Removing previous built...'
cd $HOME/.local/lib/python3.5/site-packages/
rm -r pyart
rm arm_pyart-*

# recompile
echo 'compiling'
cd $HOME/pyrad/src/pyart
python setup.py install --user