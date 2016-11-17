#!/bin/bash

# script to build pyrad and pyart
# fvj 17.11.2016

echo 'Building Pyart...'
./make_pyart.sh

echo 'Building Pyrad...'
./make_pyrad.sh
