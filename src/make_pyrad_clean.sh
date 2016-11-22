#!/bin/bash

# script to clean pyrad

# remove previous built
echo 'Removing previous built...'
cd $HOME/.local/lib/python3.5/site-packages/
rm -r pyrad
rm mch_pyrad-*

# recompile
echo 'icleaning build..'
cd $HOME/pyrad/src/pyrad_proc
python setup.py clean --all

