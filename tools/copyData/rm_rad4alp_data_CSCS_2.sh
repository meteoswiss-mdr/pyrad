#!/bin/bash

# This script removes the rad4alp data from the CSCS
# repository
# To be run in CSCS

# set permits
umask 0002

# Config
dateCmd="/bin/date"

data_destbase=/store/msrad/radar/rad4alp/rawdata/
file_type=M

while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in
        -r|--radar)
        RADAR="$2"
        OIFS=$IFS
        IFS=','
        read -r -a radar_vec <<< "$RADAR"
        IFS=$OIFS
        shift # past argument
        ;;
        -d|--day)
        DAY="$2"
        OIFS=$IFS
        IFS=','
        read -r -a day_vec <<< "$DAY"
        IFS=$OIFS
        shift # past argument
        ;;        
        -e|--res)
        res="$2"
        shift # past argument        
        ;;        
        -p|--dest_base)
        data_destbase="$2"
        shift # past argument        
        ;;    
    esac
    shift # past argument or value
done

nday=${#day_vec[@]}
nrad=${#radar_vec[@]}

for ((iday=0; iday<${nday}; iday++)); do
    yearl=$(date --date "${day_vec[${iday}]}" +"%Y")
    years=$(date --date "${day_vec[${iday}]}" +"%y")
    julday=$(date --date "${day_vec[${iday}]}" +"%j")
    
    echo "Removing day "${years}${julday}        
    
    for ((irad=0; irad<${nrad}; irad++)); do
        radar=${radar_vec[${irad}]}
        echo "Removing radar "${radar}        
        
        if [ "${res}" == "H" ]
        then            
            # remove path for polar data
            rm -rf ${data_destbase}${years}${julday}/M${res}${radar}${years}${julday}/
            rm -rf ${data_destbase}${years}${julday}/P${res}${radar}${years}${julday}/        
            
            path_ML=${data_destbase}${years}${julday}/ML${radar}${years}${julday}/
            path_PL=${data_destbase}${years}${julday}/PL${radar}${years}${julday}/
            
            # check if there is L polar data directory. If not remove also status data
            if [ ! -d "${path_ML}" ] && [ ! -d "${path_PL}" ]; then
                rm -rf ${data_destbase}${years}${julday}/ST${radar}${years}${julday}/                
            fi
        else            
            # remove path for polar data
            rm -rf ${data_destbase}${years}${julday}/M${res}${radar}${years}${julday}/
            rm -rf ${data_destbase}${years}${julday}/P${res}${radar}${years}${julday}/
                        
            # remove hydrometeor classification data
            rm -rf ${data_destbase}${years}${julday}/YM${radar}${years}${julday}/
            
            path_MH=${data_destbase}${years}${julday}/MH${radar}${years}${julday}/
            path_PH=${data_destbase}${years}${julday}/PH${radar}${years}${julday}/
            
            # check if there is L polar data directory. If not remove also status data
            if [ ! -d "${path_MH}" ] && [ ! -d "${path_PH}" ]; then
                rm -rf ${data_destbase}${years}${julday}/ST${radar}${years}${julday}/                
            fi
        fi        
    done
    
    # remove day directory if empty
    if [ ! "$(ls -A ${data_destbase}${years}${julday})" ]; then
        rm -rf ${data_destbase}${years}${julday}
    fi
done
