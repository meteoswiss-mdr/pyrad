#!/bin/bash

umask 0002;

# get solar flux data from DRAO
LOGFILE=$HOME/log/import_solar_flux.log
bash $HOME/pyrad/tools/copyData/import_solar_flux.sh >$LOGFILE 2>&1

# rad4alp dataquality monitoring
LOGFILE=$HOME/log/rad4alp_postpro_A.log
bash $HOME/pyrad/tools/processData/pyrad_rad4alp_postpro.sh A >$LOGFILE 2>&1

LOGFILE=$HOME/log/rad4alp_postpro_D.log
bash $HOME/pyrad/tools/processData/pyrad_rad4alp_postpro.sh D >$LOGFILE 2>&1

LOGFILE=$HOME/log/rad4alp_postpro_L.log
bash $HOME/pyrad/tools/processData/pyrad_rad4alp_postpro.sh L >$LOGFILE 2>&1

LOGFILE=$HOME/log/rad4alp_postpro_P.log
bash $HOME/pyrad/tools/processData/pyrad_rad4alp_postpro.sh P >$LOGFILE 2>&1

LOGFILE=$HOME/log/rad4alp_postpro_W.log
bash $HOME/pyrad/tools/processData/pyrad_rad4alp_postpro.sh W >$LOGFILE 2>&1
