#!/bin/bash

# script to build the MeteoSwiss Pyart
# fvj 17.11.2016
# sue 28.03.2018

# fetch python version in use
py=$(python --version)
pyvers=${py:7:3}

# remove previous built
echo 'Removing previous built...'
cd $HOME/.local/lib/python${pyvers}/site-packages/
rm -r pyart
rm arm_pyart-*

# recompile
echo 'compiling'
cd $HOME/pyrad/src/pyart
python setup.py install --user