#!/bin/bash
# This script controls the fetching of rad4alp Cartesian data from the
# archives and the data processing using pyrad
#
#
# Parameters:
# -d --day : days to process with format YYYYMMDD.
#            If more than one separate with a comma
# --start_time : start time of the processing each day
# --end_time : end time of the processing each day
# -c --cfgfile: config files used in the processing.
#               If more than one separate with a comma
# -p --cfgpath: path of the config files used in the processing.
#               If more than one separate with a comma
#               If more than one there must be as many as config files
# -prod: Name of the products the data of which should be retrieved.
# -i --info: Info string to be added in the name of some output products
#            If more than one config file is used separate with a comma.
#            To be used there must be as many as config files
# --get_data : if set, get product data. Default 1
# --rm_data : if set, remove product data after processing if data was fetched
#             from archive. Default 0
# --data_destbase : rad4alp data destination base path.
#                   Default /store/msrad/radar/rad4alp/rawdata/
# --MP_DSET : Parallelize dataset generation. Default 0
# --MP_PROD : Parallelize product generation. Default 0
# --MP_PROF : Profile parallelization. Default 0
umask 0002;
export PATH="/store/msrad/utils/anaconda3-pyrad/bin/:$PATH"

pyradpath="$HOME/pyrad/src/pyrad_proc/scripts/"

# default variables
START_TIME="000001"
END_TIME="240000"

GET_DATA=1
RM_DATA=0
DATA_DESTBASE=/store/msrad/radar/rad4alp/rawdata/

info_vec='None'
cfgpath_vec="$HOME/pyrad/config/processing/"

MP_DSET=0
MP_PROD=0
MP_PROF=0
while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in
        -d|--day)
        DAY="$2"
        OIFS=$IFS
        IFS=','
        read -r -a day_vec <<< "$DAY"
        IFS=$OIFS
        shift # past argument
        ;;
        --start_time)
        START_TIME="$2"
        shift # past argument
        ;;
        --end_time)
        END_TIME="$2"
        shift # past argument
        ;;
        -c|--cfgfile)
        CFGFILE="$2"
        OIFS=$IFS
        IFS=','
        read -r -a cfgfile_vec <<< "$CFGFILE"
        IFS=$OIFS
        shift # past argument
        ;;
        -p|--cfgpath)
        CFGPATH="$2"
        OIFS=$IFS
        IFS=','
        read -r -a cfgpath_vec <<< "$CFGPATH"
        IFS=$OIFS
        shift # past argument
        ;;
        --prod)
        PROD="$2"
        shift # past argument
        ;;
        -i|--info)
        INFO="$2"
        OIFS=$IFS
        IFS=','
        read -r -a info_vec <<< "$INFO"
        IFS=$OIFS
        shift # past argument
        ;;
        --get_data)
        GET_DATA="$2"
        shift # past argument
        ;;
        --rm_data)
        RM_DATA="$2"
        shift # past argument
        ;;
        --data_destbase)
        DATA_DESTBASE="$2"
        shift # past argument
        ;;
        --MP_DSET)
        MP_DSET="$2"
        shift # past argument
        ;;
        --MP_PROD)
        MP_PROD="$2"
        shift # past argument
        ;;
        --MP_PROF)
        MP_PROF="$2"
        shift # past argument
        ;;
    esac
    shift # past argument or value
done

nday=${#day_vec[@]}

ncfg=${#cfgfile_vec[@]}
nprod=${#prod_vec[@]}
ninfo=${#info_vec[@]}
ncfgpath=${#cfgpath_vec[@]}


for ((iday=0; iday<${nday}; iday++)); do
    procday_start=`date +%s`
    DAY=${day_vec[${iday}]}

    if [ "$GET_DATA" -eq 1 ]; then
        echo $(date --utc)" Getting data ${PROD} on day ${DAY}"
        proc_start=`date +%s`

        # Copy data from CSCS repository into the folder structure to process it
        LOGFILE=$HOME/log/${DAY}_get_rad4alp_cartesian_data_CSCS.log
        bash $HOME/pyrad/tools/copyData/get_rad4alp_cartesian_data_CSCS_2.sh --prod $PROD -d $DAY --start_time ${START_TIME} --end_time ${END_TIME} -p $DATA_DESTBASE >$LOGFILE 2>$LOGFILE

        proc_end=`date +%s`
        runtime=$((proc_end-proc_start))
        echo $(date --utc)" Finished getting data $PROD on day ${DAY}. Run time: ${runtime} s"
    fi

    for ((icfg=0; icfg<${ncfg}; icfg++)); do
        proccfg_start=`date +%s`

        CFGFILE=${cfgfile_vec[${icfg}]}
        IFS='.' read -r CFGFILE_BASE TERMINATION <<< "$CFGFILE"

        if [ "$ncfg" -eq "$ninfo" ];then
            INFO=${info_vec[${icfg}]}
        else
            INFO='None'
        fi
        if [ "$ncfgpath" -eq "$ncfg" ];then
            CFGPATH=${cfgpath_vec[${icfg}]}
        else
            CFGPATH="$HOME/pyrad/config/processing/"
        fi

        echo $(date --utc)" ${CFGFILE} Processing config file ${CFGFILE} on day ${DAY}"
        proc_start=`date +%s`

        # process data
        source activate pyrad

        LOGFILE=$HOME/log/${DAY}_${CFGFILE_BASE}.log
        cd ${pyradpath}

        # check start and end time
        hour_end=$(echo ${END_TIME} | cut -c1-2)
        if [ "$hour_end" -ge 24 ]; then
            day_end=$(date -d "$(date --date "${DAY}")+1 day" +"%Y%m%d")
            hour_end=$[${hour_end}-24]

            hour_end=$(printf %02d ${hour_end})
            min_end=$(echo ${END_TIME} | cut -c3-4)
            sec_end=$(echo ${END_TIME} | cut -c5-6)

            time_end=${day_end}${hour_end}${min_end}${sec_end}
        else
            time_end=${DAY}${END_TIME}
        fi
        time_start=${DAY}${START_TIME}

        python -u main_process_data.py ${CFGFILE} --starttime ${time_start} --endtime ${time_end} -i ${INFO} --cfgpath ${CFGPATH} --MULTIPROCESSING_DSET ${MP_DSET} --MULTIPROCESSING_PROD ${MP_PROD} --PROFILE_MULTIPROCESSING ${MP_PROF} >$LOGFILE 2>&1

        source deactivate

        proc_end=`date +%s`
        runtime=$((proc_end-proc_start))
        echo $(date --utc)" ${CFGFILE} Finished Processing config file ${CFGFILE} on day ${DAY}. Run time: ${runtime} s"

    done
    
    if [ "$RM_DATA" -eq 1 ] && [ "$GET_DATA" -eq 1 ]; then
        echo $(date --utc)" ${CFGFILE} Removing data ${PROD} on day ${DAY}"
        proc_start=`date +%s`
    
        LOGFILE=$HOME/log/${DAY}_rm_rad4alp_data_CSCS.log
        bash $HOME/pyrad/tools/copyData/rm_rad4alp_cartesian_data_CSCS_2.sh --prod $PROD -d $DAY --start_time ${START_TIME} --end_time ${END_TIME} -p $DATA_DESTBASE >$LOGFILE 2>$LOGFILE
    
        proc_end=`date +%s`
        runtime=$((proc_end-proc_start))
        echo $(date --utc)" ${CFGFILE} Finished removing data ${PROD} on day ${DAY}. Run time: ${runtime} s"
    fi
    
    proccfg_end=`date +%s`
    runtime=$((proccfg_end-proccfg_start))
    echo $(date --utc)" ${CFGFILE} ${DAY} Run time: ${runtime} s"
    
    
    procday_end=`date +%s`
    runtime=$((procday_end-procday_start))
    echo $(date --utc)" ${DAY} Run time: ${runtime} s"
done
