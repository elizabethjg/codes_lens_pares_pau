#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

##python -u forpair_profile.py -pcat '_zspec_best' -sample 'new' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'new_LrM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Lratio_min 0.5
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'new_Lrm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Lratio_max 0.5
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'new_zM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32  -z_min 0.4
##python -u forpair_profile.py -pcat '_zspec_best' -sample 'new_zm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32  -z_max 0.4
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'new' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32  
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'new_LrM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Lratio_min 0.5
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'new_Lrm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -Lratio_max 0.5
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'new_zM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -z_min 0.4
##python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'new_zm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -z_max 0.4
python -u forpair_profile.py -pcat '_zspec_best' -sample 'new' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'all' 
python -u forpair_profile.py -pcat '_zspec_best' -sample 'new_LrM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'all' -Lratio_min 0.5
python -u forpair_profile.py -pcat '_zspec_best' -sample 'new_Lrm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'all' -Lratio_max 0.5
python -u forpair_profile.py -pcat '_zspec_best' -sample 'new_zM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'all' -z_min 0.4
python -u forpair_profile.py -pcat '_zspec_best' -sample 'new_zm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'all' -z_max 0.4
python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'new' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'all' 
python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'new_LrM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'all' -Lratio_min 0.5
python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'new_Lrm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'all' -Lratio_max 0.5
python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'new_zM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'all' -z_min 0.4
python -u forpair_profile.py -pcat '_photo_z_2nd_run_mag_i_best' -sample 'new_zm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'all' -z_max 0.4
