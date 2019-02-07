[![Build Status](https://travis-ci.org/meteoswiss-mdr/pyrad.svg?branch=master)](https://travis-ci.org/meteoswiss-mdr/pyrad)
[![Ref doc](https://img.shields.io/badge/docs-users-4088b8.svg)](https://meteoswiss-mdr.github.io/pyrad/)

# pyrad
Python Radar Processing

# What is Pyrad?
Pyrad is a real-time data processing framework developed by MeteoSwiss. The framework is
aimed at processing and visualizing data from individual Swiss weather radars both off-line and in
real time. It is written in the Python language. The framework is version controlled and automatic
documentation is generated based on doc-strings. It is capable of ingesting data from all the
weather radars in Switzerland, namely the operational MeteoSwiss C-band rad4alp radar network,
the MeteoSwiss X-band DX50 radar and the EPFL MXPol radar and radar data in the OPERA file format.

The processing flow is controlled by 3 simple configuration files. Multiple levels of processing can
be performed. At each level new datasets (e.g. attenuation corrected reflectivity) are created which
can be stored in a file and/or used in the next processing level (for example, creating a rainfall rate
dataset from the corrected reflectivity). Multiple products can be generated from each dataset (e.g.
PPI, RHI images, histograms, etc.). In the off-line mode, data from multiple radars can be ingested
in order to obtain products such as the inter-comparison of reflectivity values at co-located range
gates.

The framework is able to ingest polarimetric and Doppler radar moments as well as auxiliary data
such as numerical weather prediction parameters (e.g. temperature, wind speed, etc.), DEM-based
visibility and data used in the generation of the products such as rain gauge measurements,
disdrometer measurements, solar flux, etc.

The signal processing and part of the data visualization is performed by a [MeteoSwiss developed version of the Py-ART radar toolkit](https://github.com/meteoswiss-mdr/pyart) which contains enhanced features. MeteoSwiss regularly contributes back to the [main Py-ART branch](https://github.com/ARM-DOE/pyart) once a new functionality has been thoroughly tested and it is considered of interest for the broad weather radar community.

The capabilities of the processing framework include various forms of echo classification and
filtering, differential phase and specific differential phase estimation, attenuation correction, data
quality monitoring, multiple rainfall rate algorithms, etc. In addition time series of data in points,
regions or trajectories of interest can be extracted and comparisons can be performed with other
sensors. This is particularly useful when performing measurement campaigns where remote
sensing retrievals are validated with in-situ airplane or ground-based measurements.

# Installation
To install Pyrad and its submodules please have a look at the Pyrad user manual in pyrad/doc/pyrad_user_manual.pdf

# Use
Before using it have a look at the cookbook in pyrad/doc/pyrad-framework-cookbook/DataProcessing.pdf

For details on the implemented functions check the pyrad library reference:

For users: pyrad/doc/pyrad_library_reference_users.pdf

For developers: pyrad/doc/pyrad_library_reference_dev.pdf


Example configuration files can be found in the directory: pyrad/config/processing/

To use Pyrad for data quality monitoring check the report: pyrad_monitoring_fvj.pdf

# Development
We welcome contributions, suggestions of developments and bug reports.

The process to contribute by partners external to MeteoSwiss is described in pyrad/doc/pyrad_user_manual.pdf

# Disclaimer
The software is still in a development stage. Please let us know if you would like to test it.

MeteoSwiss cannot be held responsible for errors in the code or problems that could arise from its use.

