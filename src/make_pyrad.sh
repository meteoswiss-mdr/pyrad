#!/bin/bash

# script to build pyrad
# fvj 17.11.2016

# remove previous built
echo 'Removing previous built...'
cd $HOME/.local/lib/python3.5/site-packages/
rm -r pyrad
rm mch_pyrad-*

# recompile
echo 'compiling...'
cd $HOME/pyrad/src/pyrad_proc
python setup.py install --user