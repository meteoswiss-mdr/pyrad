#!/bin/bash

umask 0002;

# rad4alp bird profile generation
LOGFILE=/srn/scratch/log/rad4alp_birds_A.log
bash $HOME/pyrad/tools/processData/pyrad_rad4alp_birds.sh A >$LOGFILE 2>&1

LOGFILE=/srn/scratch/log/rad4alp_birds_D.log
bash $HOME/pyrad/tools/processData/pyrad_rad4alp_birds.sh D >$LOGFILE 2>&1

LOGFILE=/srn/scratch/log/rad4alp_birds_L.log
bash $HOME/pyrad/tools/processData/pyrad_rad4alp_birds.sh L >$LOGFILE 2>&1

LOGFILE=/srn/scratch/log/rad4alp_birds_P.log
bash $HOME/pyrad/tools/processData/pyrad_rad4alp_birds.sh P >$LOGFILE 2>&1

LOGFILE=/srn/scratch/log/rad4alp_birds_W.log
bash $HOME/pyrad/tools/processData/pyrad_rad4alp_birds.sh W >$LOGFILE 2>&1
