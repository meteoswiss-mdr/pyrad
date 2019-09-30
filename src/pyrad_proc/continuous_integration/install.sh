#!/bin/bash

# set -e
# use next line to debug this script
# set -x

sudo apt-get update
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;

chmod +x miniconda.sh
./miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda
# Useful for debugging any issues with conda
conda info -a

# creation of conda environment and installation of dependencies
conda create -q -n test-environment python=$PYTHON_VERSION
source activate test-environment

# Dependencies installation:
# Py-ART required dependencies:
# - numpy
# - scipy
# - matplotlib
# - netcdf4
# Py-ART optional dependencies:
# - h5py pytest basemap cartopy gdal trmm_rsl wradlib imageio xarray
# wradlib optional dependencies:
# - xmltodict
# pyrad optional dependencies:
# - pandas shapely dask bokeh memory_profiler
conda install -c https://conda.binstar.org/jjhelmus trmm_rsl
if [[ "$PYTHON_VERSION" == "3.5" ]]; then
    # There are problems with the geos library so do not install anything that uses it
    conda install -c conda-forge numpy scipy matplotlib netcdf4 h5py pytest pandas dask bokeh memory_profiler
else
    conda install -c conda-forge numpy scipy matplotlib netcdf4 h5py pytest basemap cartopy gdal wradlib xmltodict pandas shapely dask bokeh memory_profiler imageio xarray scikit-learn
fi

# export global variables
export RSL_PATH="$HOME/miniconda/envs/test-environment"
export GDAL_DATA="$HOME/miniconda/share/gdal"
# - export PYART_CONFIG="$HOME/pyrad/config/pyart/mch_config.py"
# - export METRANETLIB_PATH=""

# installation of pyart and pyrad
cd src/pyart
python setup.py install --user
cd ../pyrad_proc
python setup.py install --user
cd
