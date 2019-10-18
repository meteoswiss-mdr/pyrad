#!/bin/bash

umask 0002;

# High resolution

## windmill 1 (ppi_nx85215) Speed greater than 0 wet
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori0.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori0_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori45.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori45_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori90.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori90_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori135.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori135_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori180.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori180_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori225.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori225_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori270.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori270_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori315.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori315_speedGT0_dBuZ_GE10
#
## windmill 1 (ppi_nx85215) Speed greater than 0 dry
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori0.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori0_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori45.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori45_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori90.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori90_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori135.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori135_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori180.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori180_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori225.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori225_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori270.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori270_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori315.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori315_speedGT0_dBuZ_LT10
#
## windmill 1 (ppi_nx85215) Speed = 0 dry
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori0.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori0_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori45.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori45_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori90.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori90_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori135.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori135_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori180.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori180_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori225.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori225_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori270.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori270_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span45.0_ori315.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori315_speedLE0_dBuZ_LT10
#
#
## windmill 2 (ppi_nx85213) Speed greater than 0 wet
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori0.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori0_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori45.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori45_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori90.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori90_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori135.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori135_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori180.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori180_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori225.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori225_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori270.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori270_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori315.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori315_speedGT0_dBuZ_GE10
#
## windmill 2 (ppi_nx85213) Speed greater than 0 dry
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori0.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori0_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori45.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori45_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori90.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori90_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori135.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori135_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori180.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori180_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori225.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori225_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori270.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori270_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori315.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori315_speedGT0_dBuZ_LT10
#
## windmill 2 (ppi_nx85213) Speed = 0 dry
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori0.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori0_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori45.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori45_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori90.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori90_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori135.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori135_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori180.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori180_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori225.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori225_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori270.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori270_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span45.0_ori315.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori315_speedLE0_dBuZ_LT10
#
#
## windmill 3 (ppi_nx85214) Speed greater than 0 wet
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori0.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori0_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori45.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori45_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori90.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori90_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori135.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori135_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori180.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori180_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori225.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori225_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori270.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori270_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori315.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori315_speedGT0_dBuZ_GE10
#
## windmill 3 (ppi_nx85214) Speed greater than 0 dry
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori0.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori0_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori45.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori45_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori90.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori90_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori135.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori135_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori180.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori180_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori225.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori225_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori270.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori270_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori315.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori315_speedGT0_dBuZ_LT10
#
#
## windmill 3 (ppi_nx85214) Speed = 0 dry
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori0.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori0_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori45.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori45_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori90.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori90_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori135.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori135_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori180.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori180_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori225.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori225_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori270.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori270_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span45.0_ori315.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori315_speedLE0_dBuZ_LT10


# RHI

## windmill 1 (rhi_nx85215) Speed greater than 0 wet
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori0.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori0_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori45.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori45_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori90.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori90_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori135.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori135_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori180.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori180_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori225.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori225_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori270.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori270_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori315.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori315_speedGT0_dBuZ_GE10
#
## windmill 1 (rhi_nx85215) Speed greater than 0 dry
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori0.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori0_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori45.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori45_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori90.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori90_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori135.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori135_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori180.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori180_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori225.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori225_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori270.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori270_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori315.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori315_speedGT0_dBuZ_LT10
#
## windmill 1 (rhi_nx85215) Speed = 0 dry
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori0.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori0_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori45.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori45_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori90.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori90_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori135.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori135_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori180.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori180_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori225.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori225_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori270.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori270_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span45.0_ori315.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori315_speedLE0_dBuZ_LT10
#
#
## windmill 2 (rhi_nx85213) Speed greater than 0 wet
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori0.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori0_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori45.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori45_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori90.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori90_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori135.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori135_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori180.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori180_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori225.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori225_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori270.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori270_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori315.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori315_speedGT0_dBuZ_GE10
#
## windmill 2 (rhi_nx85213) Speed greater than 0 dry
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori0.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori0_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori45.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori45_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori90.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori90_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori135.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori135_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori180.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori180_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori225.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori225_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori270.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori270_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori315.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori315_speedGT0_dBuZ_LT10
#
## windmill 2 (rhi_nx85213) Speed = 0 dry
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori0.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori0_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori45.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori45_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori90.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori90_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori135.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori135_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori180.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori180_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori225.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori225_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori270.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori270_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span45.0_ori315.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori315_speedLE0_dBuZ_LT10
#
## windmill 3 (rhi_nx85214) Speed greater than 0 wet
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori0.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori0_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori45.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori45_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori90.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori90_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori135.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori135_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori180.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori180_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori225.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori225_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori270.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori270_speedGT0_dBuZ_GE10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori315.0_speedGT0.0_dBuZ_GE10.0.csv --trajtype proc_periods -i span45_ori315_speedGT0_dBuZ_GE10
#
## windmill 3 (rhi_nx85214) Speed greater than 0 dry
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori0.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori0_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori45.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori45_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori90.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori90_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori135.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori135_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori180.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori180_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori225.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori225_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori270.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori270_speedGT0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori315.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori315_speedGT0_dBuZ_LT10
#
## windmill 3 (rhi_nx85214) Speed = 0 dry
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori0.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori0_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori45.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori45_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori90.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori90_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori135.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori135_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori180.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori180_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori225.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori225_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori270.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori270_speedLE0_dBuZ_LT10
#python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span45.0_ori315.0_speedLE0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span45_ori315_speedLE0_dBuZ_LT10


# HR

# windmill 1 (ppi_nx85215) Speed greater than 0 dry
python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span10.0_ori74.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span10_ori74_speedGT0_dBuZ_LT10
python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM1.txt -t /users/jfigui/windmills_params/ppi_nx85215_span10.0_ori251.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span10_ori251_speedGT0_dBuZ_LT10

# windmill 2 (ppi_nx85213) Speed greater than 0 dry
python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span10.0_ori55.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span10_ori55_speedGT0_dBuZ_LT10
python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM2.txt -t /users/jfigui/windmills_params/ppi_nx85213_span10.0_ori214.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span10_ori214_speedGT0_dBuZ_LT10

# windmill 3 (ppi_nx85214) Speed greater than 0 dry
python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span10.0_ori89.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span10_ori89_speedGT0_dBuZ_LT10
python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_hr_WM3.txt -t /users/jfigui/windmills_params/ppi_nx85214_span10.0_ori250.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span10_ori250_speedGT0_dBuZ_LT10


# RHI

# windmill 1 (rhi_nx85215) Speed greater than 0 dry
python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span10.0_ori78.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span10_ori78_speedGT0_dBuZ_LT10
python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM1.txt -t /users/jfigui/windmills_params/rhi_nx85215_span10.0_ori246.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span10_ori246_speedGT0_dBuZ_LT10

# windmill 2 (rhi_nx85213) Speed greater than 0 dry
python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span10.0_ori54.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span10_ori54_speedGT0_dBuZ_LT10
python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM2.txt -t /users/jfigui/windmills_params/rhi_nx85213_span10.0_ori213.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span10_ori213_speedGT0_dBuZ_LT10

# windmill 3 (rhi_nx85214) Speed greater than 0 dry
python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span10.0_ori89.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span10_ori89_speedGT0_dBuZ_LT10
python -u /users/jfigui/pyrad/src/pyrad_proc/scripts/main_process_data.py cscs_mals_sha_windmills_rhi_WM3.txt -t /users/jfigui/windmills_params/rhi_nx85214_span10.0_ori250.0_speedGT0.0_dBuZ_LT10.0.csv --trajtype proc_periods -i span10_ori250_speedGT0_dBuZ_LT10
