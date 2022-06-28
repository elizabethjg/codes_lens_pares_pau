#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

##python -u forpair_profile.py -pcat '_zspec_best' -sample 'new' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'new_LrM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Lratio_min 0.8
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'new_Lrm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Lratio_max 0.8
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'new_zM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32  -z_min 0.4
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'new_zm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32  -z_max 0.4
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'new' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32  
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'new_LrM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Lratio_min 0.8
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'new_Lrm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Lratio_max 0.8
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'new_zM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -z_min 0.4
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'new_zm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -z_max 0.4
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'mh' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' 
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'mh_LrM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Lratio_min 0.8
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'mh_Lrm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Lratio_max 0.8
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'mh_zM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -z_min 0.4
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'mh_zm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -z_max 0.4
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' 
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_LrM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Lratio_min 0.8
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_Lrm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Lratio_max 0.8
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_zM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -z_min 0.4
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_zm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -z_max 0.4
##python -u forpair_profile.py -pcat '_zspec' -sample 'mh' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' 
##python -u forpair_profile.py -pcat '_zspec' -sample 'mh_LrM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Lratio_min 0.8
##python -u forpair_profile.py -pcat '_zspec' -sample 'mh_Lrm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Lratio_max 0.8
##python -u forpair_profile.py -pcat '_zspec' -sample 'mh_zM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -z_min 0.4
##python -u forpair_profile.py -pcat '_zspec' -sample 'mh_zm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -z_max 0.4
python -u forpair_profile.py -pcat '_vane' -sample 'mh' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' 
python -u forpair_profile.py -pcat '_vane' -sample 'mh_LrM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Lratio_min 0.8
python -u forpair_profile.py -pcat '_vane' -sample 'mh_Lrm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Lratio_max 0.8
python -u forpair_profile.py -pcat '_vane' -sample 'mh_zM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -z_min 0.4
python -u forpair_profile.py -pcat '_vane' -sample 'mh_zm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -z_max 0.4
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'mh' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' 
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'mh_LrM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Lratio_min 0.8
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'mh_Lrm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Lratio_max 0.8
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'mh_zM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -z_min 0.4
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'mh_zm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -z_max 0.4
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' 
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_LrM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Lratio_min 0.8
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_Lrm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Lratio_max 0.8
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_zM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -z_min 0.4
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_zm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -z_max 0.4
##python -u forpair_profile.py -pcat '_zspec' -sample 'mh' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' 
##python -u forpair_profile.py -pcat '_zspec' -sample 'mh_LrM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Lratio_min 0.8
##python -u forpair_profile.py -pcat '_zspec' -sample 'mh_Lrm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Lratio_max 0.8
##python -u forpair_profile.py -pcat '_zspec' -sample 'mh_zM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -z_min 0.4
##python -u forpair_profile.py -pcat '_zspec' -sample 'mh_zm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -z_max 0.4
python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' 
python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_LrM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Lratio_min 0.8
python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_Lrm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Lratio_max 0.8
python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_zM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -z_min 0.4
python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_zm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -z_max 0.4
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'mh_red' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -colorcut 'red'
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'mh_blue' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -colorcut 'blue'
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'mh_MM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Mtot_min -21.0
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'mh_Mm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Mtot_max -21.0
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_red' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -colorcut 'red'
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_blue' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -colorcut 'blue'
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_MM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Mtot_min -21.0
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_Mm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Mtot_max -21.0
##python -u forpair_profile.py -pcat '_zspec' -sample 'mh_red' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -colorcut 'red'
##python -u forpair_profile.py -pcat '_zspec' -sample 'mh_blue' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -colorcut 'blue'
##python -u forpair_profile.py -pcat '_zspec' -sample 'mh_MM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Mtot_min -21.0
##python -u forpair_profile.py -pcat '_zspec' -sample 'mh_Mm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Mtot_max -21.0
python -u forpair_profile.py -pcat '_vane' -sample 'mh_red' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -colorcut 'red'
python -u forpair_profile.py -pcat '_vane' -sample 'mh_blue' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -colorcut 'blue'
python -u forpair_profile.py -pcat '_vane' -sample 'mh_MM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Mtot_min -21.0
python -u forpair_profile.py -pcat '_vane' -sample 'mh_Mm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Mtot_max -21.0
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'mh_red' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -colorcut 'red'
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'mh_blue' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -colorcut 'blue'
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'mh_MM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Mtot_min -21.0
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'mh_Mm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Mtot_max -21.0
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_red' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -colorcut 'red'
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_blue' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -colorcut 'blue'
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_MM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Mtot_min -21.0
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_Mm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Mtot_max -21.0
##python -u forpair_profile.py -pcat '_zspec' -sample 'mh_red' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -colorcut 'red'
##python -u forpair_profile.py -pcat '_zspec' -sample 'mh_blue' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -colorcut 'blue'
##python -u forpair_profile.py -pcat '_zspec' -sample 'mh_MM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Mtot_min -21.0
##python -u forpair_profile.py -pcat '_zspec' -sample 'mh_Mm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Mtot_max -21.0
python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_red' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -colorcut 'red'
python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_blue' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -colorcut 'blue'
python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_MM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Mtot_min -21.0
python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_Mm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1' -Mtot_max -21.0
