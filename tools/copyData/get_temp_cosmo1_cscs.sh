#!/bin/bash
# This script: 
# Gets COSMO-1 temperature file from CSCS repository and puts it in the right
# folder to be processed in CSCS 
# created 23.03.2017 (fvj)

# set permits
umask 0002

# Config
dateCmd="/bin/date"
cosmobaseraw='/store/msrad/cosmo/cosmo1/TEMP/raw1/'
cosmopathcscs='/store/s83/owm/COSMO-1/ORDERS/MDR/'

#module load idl
. /apps/albis/itt/idl/idl84/inst/idl84/bin/idl_setup.bash

export IDL_OTLDIR=/store/msrad/idl
export IDL_STARTUP=$IDL_OTLDIR/im_setup.idl;export IDL_STARTUP

export IDL_DIR=/apps/albis/itt/idl/idl84/inst/idl84/
export LD_LIBRARY_PATH=$HOME/malsgit/src/libDX50/lib:$LD_LIBRARY_PATH

export IDL_NO_XWINDOWS=1
export DISPLAY=

execpath=$HOME/malsgit/src/idl/cosmo/
idlpath=/apps/albis/itt/idl/idl84/inst/idl84/bin/idl

# defaults
radar="None"
res="L"
while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in        
        -d|--day)
        DAY="$2"
        OIFS=$IFS
        IFS=','
        read -r -a date_vec <<< "$DAY"
        IFS=$OIFS
        shift # past argument
        ;;
        -t|--time)
        RUN="$2"
        OIFS=$IFS
        IFS=','
        read -r -a hour_run <<< "$RUN"
        IFS=$OIFS
        shift # past argument
        ;;        
        -r|--radar)
        radar="$2"
        shift # past argument
        ;;
        -e|--res)
        res="$2"
        shift # past argument
        ;;
    esac
    shift # past argument or value
done

nday=${#date_vec[@]}
nhour=${#hour_run[@]}

# Log
date

for ((iday=0; iday<${nday}; iday++)); do
    datedir=$(${dateCmd} --date "${date_vec[${iday}]}" +"%Y-%m-%d")

    for ((ihour=0; ihour<${nhour}; ihour++)); do		
        # import COSMO temperature data from cscs
        cosmoFileRaw=cosmo-1_MDR_3D_${date_vec[${iday}]}${hour_run[${ihour}]}.nc
        cd ${cosmobaseraw}
        
        echo "Importing COSMO file "${cosmoFileRaw}" from cscs repository"
		
        mkdir -p ${datedir}
        cp ${cosmopathcscs}${cosmoFileRaw} ${datedir}
        
        if [ $radar != "None" ]
        then                    
            years=$(date --date "${date_vec[${iday}]}" +"%y")
            julday=$(date --date "${date_vec[${iday}]}" +"%j")
                    
            temp_cosmo_radar=${radar}
            temp_cosmo_res=${res}
            temp_cosmo_year=$(date --date "${date_vec[${iday}]}" +"%y")
            temp_cosmo_day=$(date --date "${date_vec[${iday}]}" +"%j")
            temp_cosmo_hour=${hour_run[${ihour}]}

            export temp_cosmo_radar
            export temp_cosmo_res
            export temp_cosmo_year
            export temp_cosmo_day
            export temp_cosmo_hour

            cd ${execpath}
            ${idlpath} -rt=temp_cosmo1_cscs.run
        fi
    done
done

# Log
echo "All done!"
date
