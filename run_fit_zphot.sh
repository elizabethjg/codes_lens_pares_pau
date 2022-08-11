#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

##python -u fit_pairs_profile.py -ncores 5 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 15 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_LrM2_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 15 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_Lrm2_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 5 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_zm_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 5 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_zM_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 4 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_MM_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 4 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_Mm_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 4 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_red_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 4 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_blue_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 4 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_M1_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 4 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_M2_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 4 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_M3_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 4 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_M4_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
wait
