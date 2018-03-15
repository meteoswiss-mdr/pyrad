#!/bin/bash

umask 0002;

# rad4alp ground clutter monitoring
LOGFILE=/srn/scratch/log/rad4alp_gc_A.log
bash $HOME/pyrad/tools/processData/pyrad_rad4alp_gc_monitoring.sh A >$LOGFILE 2>&1

LOGFILE=/srn/scratch/log/rad4alp_gc_D.log
bash $HOME/pyrad/tools/processData/pyrad_rad4alp_gc_monitoring.sh D >$LOGFILE 2>&1

LOGFILE=/srn/scratch/log/rad4alp_gc_L.log
bash $HOME/pyrad/tools/processData/pyrad_rad4alp_gc_monitoring.sh L >$LOGFILE 2>&1

LOGFILE=/srn/scratch/log/rad4alp_gc_P.log
bash $HOME/pyrad/tools/processData/pyrad_rad4alp_gc_monitoring.sh P >$LOGFILE 2>&1

LOGFILE=/srn/scratch/log/rad4alp_gc_W.log
bash $HOME/pyrad/tools/processData/pyrad_rad4alp_gc_monitoring.sh W >$LOGFILE 2>&1
