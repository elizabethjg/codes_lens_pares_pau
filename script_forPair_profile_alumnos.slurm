#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32  
python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_LrM2' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Lratio_min 0.5
python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_Lrm2' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Lratio_max 0.5
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_zM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -z_min 0.4
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_zm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -z_max 0.4
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_red' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -colorcut 'red'
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_blue' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -colorcut 'blue'
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_MM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Mtot_min -21.0
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_Mm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Mtot_max -21.0
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_M1' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Mtot_max -21.09
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_M2' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Mtot_min -21.09 -Mtot_max -20.65
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_M3' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Mtot_min -20.65 -Mtot_max -20.28
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_M4' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Mtot_min -20.28
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_M_1' -nbins 15 -RIN 100 -ROUT 10000 -ncores 64 -Mtot_max -20.93
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_M_2' -nbins 15 -RIN 100 -ROUT 10000 -ncores 64 -Mtot_min -20.93 -Mtot_max -20.40
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i' -sample 'mh_M_3' -nbins 15 -RIN 100 -ROUT 10000 -ncores 64 -Mtot_min -20.40
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32  
python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_LrM2' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Lratio_min 0.5
python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_Lrm2' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Lratio_max 0.5
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_zM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -z_min 0.4
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_zm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -z_max 0.4
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_red' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -colorcut 'red'
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_blue' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -colorcut 'blue'
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_MM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Mtot_min -21.0
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_Mm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Mtot_max -21.0
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_M1' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Mtot_max -21.32
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_M2' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Mtot_min -21.32 -Mtot_max -20.91
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_M3' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Mtot_min -20.91 -Mtot_max -20.41
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_M4' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Mtot_min -20.41
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_M_1' -nbins 15 -RIN 100 -ROUT 10000 -ncores 64 -Mtot_max -21.17
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_M_2' -nbins 15 -RIN 100 -ROUT 10000 -ncores 64 -Mtot_min -21.17 -Mtot_max -20.59
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'mh_M_3' -nbins 15 -RIN 100 -ROUT 10000 -ncores 64 -Mtot_min -20.59 

