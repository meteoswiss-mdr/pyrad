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
    
    yearl_end=$(date -d "$(date --date "${day_vec[${iday}]}")+1 day" +"%Y")
    years_end=$(date -d "$(date --date "${day_vec[${iday}]}")+1 day" +"%y")
    julday_end=$(date -d "$(date --date "${day_vec[${iday}]}")+1 day" +"%j")
    
    echo "Removing day "${years}${julday}        
    
    for ((irad=0; irad<${nrad}; irad++)); do
        radar=${radar_vec[${irad}]}
        echo "Removing radar "${radar}        
        
        if [ "${res}" == "H" ]
        then            
            # remove polar data from 00:05 to 23:55
            find ${data_destbase}${years}${julday}/M${res}${radar}${years}${julday}/ -type f -not -name M${res}${radar}${years}${julday}0000*.* -delete
            find ${data_destbase}${years}${julday}/P${res}${radar}${years}${julday}/ -type f -not -name P${res}${radar}${years}${julday}0000*.* -delete
            
            # remove directory if empty
            if [ ! "$(ls -A ${data_destbase}${years}${julday}/M${res}${radar}${years}${julday})" ]; then
                rm -rf ${data_destbase}${years}${julday}/M${res}${radar}${years}${julday}
            fi
            if [ ! "$(ls -A ${data_destbase}${years}${julday}/P${res}${radar}${years}${julday})" ]; then
                rm -rf ${data_destbase}${years}${julday}/P${res}${radar}${years}${julday}
            fi
                        
            # remove polar data at 00:00 the next day
            rm -f ${data_destbase}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end}0000*.*
            rm -f ${data_destbase}${years_end}${julday_end}/P${res}${radar}${years_end}${julday_end}/P${res}${radar}${years_end}${julday_end}0000*.*
            
            # remove directory if empty
            if [ ! "$(ls -A ${data_destbase}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end})" ]; then
                rm -rf ${data_destbase}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end}
            fi
            if [ ! "$(ls -A ${data_destbase}${years_end}${julday_end}/P${res}${radar}${years_end}${julday_end})" ]; then
                rm -rf ${data_destbase}${years_end}${julday_end}/P${res}${radar}${years_end}${julday_end}
            fi
                        
            # check if there is L polar data directory. If not remove also status data
            path_ML=${data_destbase}${years}${julday}/ML${radar}${years}${julday}/
            path_PL=${data_destbase}${years}${julday}/PL${radar}${years}${julday}/            
            if [ ! -d "${path_ML}" ] && [ ! -d "${path_PL}" ]; then
                find ${data_destbase}${years}${julday}/ST${radar}${years}${julday}/ -type f -not -name ST${radar}${years}${julday}0000*.xml -delete                                
            fi
            # remove directory if empty
            if [ ! "$(ls -A ${data_destbase}${years}${julday}/ST${radar}${years}${julday})" ]; then
                rm -rf ${data_destbase}${years}${julday}/ST${radar}${years}${julday}
            fi
            
            # check if there is L polar data directory next day. If not remove also status data
            path_ML=${data_destbase}${years_end}${julday_end}/ML${radar}${years_end}${julday_end}/
            path_PL=${data_destbase}${years_end}${julday_end}/PL${radar}${years_end}${julday_end}/            
            if [ ! -d "${path_ML}" ] && [ ! -d "${path_PL}" ]; then
                rm -f ${data_destbase}${years_end}${julday_end}/ST${radar}${years_end}${julday_end}/ST${radar}${years_end}${julday_end}0000*.xml
            fi
            # remove directory if empty
            if [ ! "$(ls -A ${data_destbase}${years_end}${julday_end}/ST${radar}${years_end}${julday_end})" ]; then
                rm -rf ${data_destbase}${years_end}${julday_end}/ST${radar}${years_end}${julday_end}
            fi
        else            
            # remove polar data from 00:05 to 23:55
            find ${data_destbase}${years}${julday}/M${res}${radar}${years}${julday}/ -type f -not -name M${res}${radar}${years}${julday}0000*.* -delete
            find ${data_destbase}${years}${julday}/P${res}${radar}${years}${julday}/ -type f -not -name P${res}${radar}${years}${julday}0000*.* -delete
            
            # remove directory if empty
            if [ ! "$(ls -A ${data_destbase}${years}${julday}/M${res}${radar}${years}${julday})" ]; then
                rm -rf ${data_destbase}${years}${julday}/M${res}${radar}${years}${julday}
            fi
            if [ ! "$(ls -A ${data_destbase}${years}${julday}/P${res}${radar}${years}${julday})" ]; then
                rm -rf ${data_destbase}${years}${julday}/P${res}${radar}${years}${julday}
            fi
            
            # remove polar data at 00:00 the next day
            rm -f ${data_destbase}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end}0000*.*
            rm -f ${data_destbase}${years_end}${julday_end}/P${res}${radar}${years_end}${julday_end}/P${res}${radar}${years_end}${julday_end}0000*.*
            
            # remove directory if empty
            if [ ! "$(ls -A ${data_destbase}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end})" ]; then
                rm -rf ${data_destbase}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end}
            fi
            if [ ! "$(ls -A ${data_destbase}${years_end}${julday_end}/P${res}${radar}${years_end}${julday_end})" ]; then
                rm -rf ${data_destbase}${years_end}${julday_end}/P${res}${radar}${years_end}${julday_end}
            fi
                        
            # remove hydrometeor classification data
            find ${data_destbase}${years}${julday}/YM${radar}${years}${julday}/ -type f -not -name YM${radar}${years}${julday}0000*.* -delete            
            rm -f ${data_destbase}${years_end}${julday_end}/YM${radar}${years_end}${julday_end}/YM${radar}${years_end}${julday_end}0000*.*
            
            # remove directory if empty
            if [ ! "$(ls -A ${data_destbase}${years}${julday}/YM${radar}${years}${julday})" ]; then
                rm -rf ${data_destbase}${years}${julday}/YM${radar}${years}${julday}
            fi
            if [ ! "$(ls -A ${data_destbase}${years_end}${julday_end}/YM${radar}${years_end}${julday_end})" ]; then
                rm -rf ${data_destbase}${years_end}${julday_end}/YM${radar}${years_end}${julday_end}
            fi            
            
            path_MH=${data_destbase}${years}${julday}/MH${radar}${years}${julday}/
            path_PH=${data_destbase}${years}${julday}/PH${radar}${years}${julday}/
            
            # check if there is L polar data directory. If not remove also status data
            if [ ! -d "${path_MH}" ] && [ ! -d "${path_PH}" ]; then
                find ${data_destbase}${years}${julday}/ST${radar}${years}${julday}/ -type f -not -name ST${radar}${years}${julday}0000*.xml -delete
            fi
            # remove directory if empty
            if [ ! "$(ls -A ${data_destbase}${years}${julday}/ST${radar}${years}${julday})" ]; then
                rm -rf ${data_destbase}${years}${julday}/ST${radar}${years}${julday}
            fi
            
            # check if there is L polar data directory next day. If not remove also status data
            path_ML=${data_destbase}${years_end}${julday_end}/ML${radar}${years_end}${julday_end}/
            path_PL=${data_destbase}${years_end}${julday_end}/PL${radar}${years_end}${julday_end}/            
            if [ ! -d "${path_ML}" ] && [ ! -d "${path_PL}" ]; then
                rm -f ${data_destbase}${years_end}${julday_end}/ST${res}${radar}${years_end}${julday_end}/ST${res}${radar}${years_end}${julday_end}0000*.xml
            fi
            # remove directory if empty
            if [ ! "$(ls -A ${data_destbase}${years_end}${julday_end}/ST${radar}${years_end}${julday_end})" ]; then
                rm -rf ${data_destbase}${years_end}${julday_end}/ST${radar}${years_end}${julday_end}
            fi
        fi        
    done
    
    # remove day directory if empty
    if [ ! "$(ls -A ${data_destbase}${years}${julday})" ]; then
        rm -rf ${data_destbase}${years}${julday}
    fi
    if [ ! "$(ls -A ${data_destbase}${years_end}${julday_end})" ]; then
        rm -rf ${data_destbase}${years_end}${julday_end}
    fi
done
