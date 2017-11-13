#!/bin/bash

# This script gets the rad4alp data from the CSCS
# repository and puts it in the right folder.
# To be run in CSCS

# set permits
umask 0002

# Config
dateCmd="/bin/date"

phdata_origbase=/store/msrad/radar/polarHR/data/
rawdata_origbase=/store/msrad/radar/swiss/data/
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
    esac
    shift # past argument or value
done

nday=${#day_vec[@]}
nrad=${#radar_vec[@]}

for ((irad=0; irad<${nrad}; irad++)); do
    radar=${radar_vec[${irad}]}
    echo "Processing radar "${radar}
    
    for ((iday=0; iday<${nday}; iday++)); do
        yearl=$(date --date "${day_vec[${iday}]}" +"%Y")
        years=$(date --date "${day_vec[${iday}]}" +"%y")
        julday=$(date --date "${day_vec[${iday}]}" +"%j")
    
        echo "Processing day "${years}${julday}        
        # transfer polar data from CSCS to destination folder, unzip it and remove zip file
        if [ "${res}" == "H" ]
        then
            # check type of file that exists in the repository
            if [ `ls ${phdata_origbase}${yearl}/${years}${julday}/M${res}${radar}${years}${julday}.zip` ]
            then
                file_type=M
            elif [ `ls ${phdata_origbase}${yearl}/${years}${julday}/P${res}${radar}${years}${julday}.zip` ]
            then
                file_type=P
            else
                echo "No file found in repository"
                continue
            fi            
            
            # create destination path for polar data
            mkdir -p ${data_destbase}${years}${julday}/${file_type}${res}${radar}${years}${julday}/
            cd ${data_destbase}${years}${julday}/${file_type}${res}${radar}${years}${julday}/
            
            cp ${phdata_origbase}${yearl}/${years}${julday}/${file_type}${res}${radar}${years}${julday}.zip .
            unzip -o ${file_type}${res}${radar}${years}${julday}.zip
            chmod -R gu+rw ${data_destbase}${years}${julday}/${file_type}${res}${radar}${years}${julday}/
            rm ${file_type}${res}${radar}${years}${julday}.zip            
        else
            # check type of file that exists in the repository
            if [ `ls ${rawdata_origbase}${yearl}/${years}${julday}/M${res}${radar}${years}${julday}.zip` ]
            then
                file_type=M
            elif [ `ls ${rawdata_origbase}${yearl}/${years}${julday}/P${res}${radar}${years}${julday}.zip` ]
            then
                file_type=P
            else
                echo "No file found in repository"
                continue
            fi
            
            # create destination path for polar data
            mkdir -p ${data_destbase}${years}${julday}/${file_type}${res}${radar}${years}${julday}/
            cd ${data_destbase}${years}${julday}/${file_type}${res}${radar}${years}${julday}/
            
            cp ${rawdata_origbase}${yearl}/${years}${julday}/${file_type}${res}${radar}${years}${julday}.zip .
            unzip -o ${file_type}${res}${radar}${years}${julday}.zip
            chmod -R gu+rw ${data_destbase}${years}${julday}/${file_type}${res}${radar}${years}${julday}/
            rm ${file_type}${res}${radar}${years}${julday}.zip
            
            # create destination path for hydrometeor classification data
            mkdir -p ${data_destbase}${years}${julday}/YM${radar}${years}${julday}/
            cd ${data_destbase}${years}${julday}/YM${radar}${years}${julday}/
            
            cp ${phdata_origbase}${yearl}/${years}${julday}/YM${radar}${years}${julday}.zip .
            unzip -o YM${radar}${years}${julday}.zip
            chmod -R gu+rw ${data_destbase}${years}${julday}/YM${radar}${years}${julday}/
            rm YM${radar}${years}${julday}.zip
        fi
    
        # create destination path for status data
        mkdir -p ${data_destbase}${years}${julday}/ST${radar}${years}${julday}/
        cd ${data_destbase}${years}${julday}/ST${radar}${years}${julday}/
        
        # transfer status data from CSCS to destination folder, unzip it and remove zip file
        cp ${rawdata_origbase}${yearl}/${years}${julday}/ST${radar}${years}${julday}.zip .
        unzip -o ST${radar}${years}${julday}.zip
        chmod -R gu+rw ${data_destbase}${years}${julday}/ST${radar}${years}${julday}/
        rm ST${radar}${years}${julday}.zip
    done
done
