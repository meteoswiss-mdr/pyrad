#!/bin/bash

# Call the realtime processing application with arguments.
#  $1 : config file
#  $2 : log file
#  $3 : Set to 1, if the log file must be renamed.
#  $4 : Appendix to log file if $3 is true.
function postprocessing {
    if [ $4 != 0 ]; then
        # Rename logfiles (day changed)
        if [ -f $3 ]; then
            mv $3 $3$5
        fi
    fi

    # Run postproc processing
    cd ${pyradpath}
    python main_process_data.py $1 -i $2 >>$3 2>>$3
}

# Call the realtime processing application with arguments.
#  $1 : config file
#  $2 : start date of the processing
#  $3 : End date of the processing
#  $4 : log file
#  $5 : Set to 1, if the log file must be renamed.
#  $6 : Appendix to log file if $5 is true.
function dataquality {
    if [ $6 != 0 ]; then
        # Rename logfiles (day changed)
        if [ -f $5 ]; then
            mv $5 $5$7
        fi
    fi
    
    # Run postproc processing
    cd ${pyradpath}
    python main_process_data_period.py $1 $2 $3 -i $4 >>$5 2>>$5
}

# set permits
umask 0002

# input variables
RADAR=$1

echo "PROCESSING RADAR "${RADAR}

# activate pyrad environment
source /srn/analysis/anaconda3/bin/activate pyrad

proc_start=`date +%s`

pyradpath="$HOME/pyrad/src/pyrad_proc/scripts/"

# File where to save day of last cron run
POSTPROC_LASTSTATE="$HOME/postproc_pyrad/${RADAR}_laststate.txt"

# Check if new day: if yes rename logfiles
RENAME_LOGFILES=0
LOG_APPENDIX=''
TODAY=`date --utc +"%Y%m%d"`
if [ -f $POSTPROC_LASTSTATE ]; then
    LASTSTATE=`cat $POSTPROC_LASTSTATE`
    if [ $TODAY != $LASTSTATE ]; then
        RENAME_LOGFILES=1
        LOG_APPENDIX="_$LASTSTATE"
        echo $TODAY > $POSTPROC_LASTSTATE
    fi
    START_TIME=$LASTSTATE
else
    echo $TODAY > $POSTPROC_LASTSTATE
    START_TIME=$(date --date ${TODAY}'-24 hours' +"%Y%m%d")
fi
END_TIME=$(date --date ${TODAY}'-24 hours' +"%Y%m%d")

# PL Data quality
CONFIGFILE=rad4alp_dataquality_PL${RADAR}.txt
LOGFILE=$HOME/log/rad4alp_dataquality_PL${RADAR}.log
dataquality $CONFIGFILE  $START_TIME $END_TIME $RADAR $LOGFILE $RENAME_LOGFILES $LOG_APPENDIX

# # PH Data quality (sun monitoring)
# CONFIGFILE=rad4alp_dataquality_PH${RADAR}.txt
# LOGFILE=$HOME/log/rad4alp_dataquality_PH${RADAR}.log
# dataquality $CONFIGFILE  $START_TIME $END_TIME $RADAR $LOGFILE $RENAME_LOGFILES $LOG_APPENDIX

source /srn/analysis/anaconda3/bin/deactivate

proc_end=`date +%s`
runtime=$((proc_end-proc_start))

echo "Run time: ${runtime} s"
