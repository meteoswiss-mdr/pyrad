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
    python -u main_process_data.py $1 -i $2 >>$3 2>>$3
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
    python -u main_process_data_birds.py $1 $2 $3 --storepath $4 >>$5 2>>$5
}

# set permits
umask 0002

# input variables
RADAR=$1

CURRENT_TIME=$(date --utc)

# activate pyrad environment
source /srn/analysis/anaconda3/bin/activate pyrad

proc_start=`date +%s`

pyradpath="$HOME/pyrad/src/pyrad_proc/scripts/"

# File where to save day of last cron run
POSTPROC_LASTSTATE="$HOME/birds_pyrad/${RADAR}_laststate.txt"

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
END_TIME=$(date --date ${TODAY} +"%Y%m%d")

# log data
echo "PROCESSING RADAR "${RADAR}
echo "PROCESSING START TIME: "${CURRENT_TIME}
echo "START TIME OF DATA TO BE PROCESSED "${START_TIME}000001
echo "END TIME OF DATA TO BE PROCESSED "${END_TIME}000000

# PH bird detection
CONFIGFILE=rad4alp_birds_PH${RADAR}.txt
LOGFILE=/srn/scratch/log/rad4alp_pyrad_birds_PH${RADAR}.log
STOREPATH=/srn/analysis/pyrad_products/rad4alp_birds_PH${RADAR}/
dataquality $CONFIGFILE  ${START_TIME}000001 ${END_TIME}000000 $STOREPATH $LOGFILE $RENAME_LOGFILES $LOG_APPENDIX

# Transfer data to Vogelwarte
cd $STOREPATH
zipfile=r0082z.pyradbird_${RADAR}_${START_TIME}.zip
zip -r ${zipfile} $(date -d ${START_TIME} +%Y-%m-%d)
/tools/ext/wman/custom/bin/socupload -s zuegts34 -f ${zipfile}
rm ${zipfile}

source /srn/analysis/anaconda3/bin/deactivate

proc_end=`date +%s`
runtime=$((proc_end-proc_start))

CURRENT_TIME=$(date --utc)
echo "PROCESSING END TIME: "${CURRENT_TIME}
echo "Total run time: ${runtime} s"
