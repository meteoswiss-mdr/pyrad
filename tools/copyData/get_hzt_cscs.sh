#!/bin/bash

# This script gets the HZT data from the CSCS
# repository and puts it in the right folder.
# To be run in CSCS

# set permits
umask 0002

# Config
dateCmd="/bin/date"

rawdata_origbase=/store/msrad/radar/swiss/data/
data_destbase=/store/msrad/cosmo/cosmo1/HZT/

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
    esac
    shift # past argument or value
done

nday=${#day_vec[@]}
nrad=${#radar_vec[@]}

for ((iday=0; iday<${nday}; iday++)); do
    yearl=$(date --date "${day_vec[${iday}]}" +"%Y")
    years=$(date --date "${day_vec[${iday}]}" +"%y")
    julday=$(date --date "${day_vec[${iday}]}" +"%j")

    echo "Processing day "${years}${julday}

    # create destination path for polar data
    mkdir -p ${data_destbase}${years}${julday}/
    cd ${data_destbase}${years}${julday}/

    # transfer HZT from CSCS to destination folder, unzip it and remove zip file
    cp ${rawdata_origbase}${yearl}/${years}${julday}/HZT${years}${julday}.zip .
    unzip -o HZT${years}${julday}.zip
    chmod -R gu+rw ${data_destbase}${years}${julday}/
    rm HZT${years}${julday}.zip
done
