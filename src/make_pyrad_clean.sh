#!/bin/bash

# script to clean pyrad

# fetch python version in use
py=$(python --version)
pyvers=${py:7:3}

# remove previous built
echo 'Removing previous built...'
cd $HOME/.local/lib/python${pyvers}/site-packages/
rm -r pyrad
rm mch_pyrad-*

# recompile
echo 'icleaning build..'
cd $HOME/pyrad/src/pyrad_proc
python setup.py clean --all

