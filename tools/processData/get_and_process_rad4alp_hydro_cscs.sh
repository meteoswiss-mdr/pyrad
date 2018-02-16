#!/bin/bash

umask 0002;
export PATH="/store/msrad/utils/anaconda3/bin:$PATH"

pyradpath="$HOME/pyrad/src/pyrad_proc/scripts/"

cfgfile=cscs_rad4alp_hydro_PLA.txt

DAY=20170818
START_TIME=180000
END_TIME=190000

RADAR=A
RES=L
RUN_TIME=18


# Copy data from CSCS repository into the folder structure to process it
LOGFILE=$HOME/log/${DAY}_${RADAR}${RES}_get_rad4alp_data_CSCS.log
bash $HOME/pyrad/tools/copyData/get_rad4alp_data_CSCS_2.sh -r $RADAR -e $RES -d $DAY >$LOGFILE 2>$LOGFILE

# Get COSMO temperature information files
LOGFILE=$HOME/log/${DAY}_${RADAR}${RES}_temp_cosmo1_rad4alp.log
bash $HOME/pyrad/tools/copyData/get_temp_cosmo1_cscs.sh -d $DAY -t ${RUN_TIME} >$LOGFILE 2>$LOGFILE

# copy HZT files to right folder
LOGFILE=$HOME/log/${DAY}_${RADAR}${RES}_hzt_rad4alp.log
bash $HOME/pyrad/tools/copyData/get_hzt_cscs.sh -d $DAY >$LOGFILE 2>$LOGFILE

# process trajectory
source activate pyrad

LOGFILE=$HOME/log/${DAY}_${RADAR}${RES}_rad4alp_hydro.log
cd ${pyradpath}
python main_process_data.py ${cfgfile} --starttime ${DAY}${START_TIME} --endtime ${DAY}${END_TIME} >$LOGFILE 2>&1

source deactivate