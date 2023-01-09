#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

##python -u fit_pairs_profile.py -ncores 7 -folder ~/PARES-PAU/profiles_pk/ -file 'profile_pk_w1_Vane.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 7 -folder ~/PARES-PAU/profiles_pk/ -file 'profile_pk_LrM2_w1_Vane.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 7 -folder ~/PARES-PAU/profiles_pk/ -file 'profile_pk_Lrm2_w1_Vane.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 7 -folder ~/PARES-PAU/profiles_pk/ -file 'profile_pk_zm_w1_Vane.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 7 -folder ~/PARES-PAU/profiles_pk/ -file 'profile_pk_zM_w1_Vane.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 7 -folder ~/PARES-PAU/profiles_pk/ -file 'profile_pk_MM_w1_Vane.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 7 -folder ~/PARES-PAU/profiles_pk/ -file 'profile_pk_Mm_w1_Vane.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 7 -folder ~/PARES-PAU/profiles_pk/ -file 'profile_pk_red_w1_Vane.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 7 -folder ~/PARES-PAU/profiles_pk/ -file 'profile_pk_blue_w1_Vane.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 10 -folder ~/PARES-PAU/profiles_pk/ -file 'profile_pk_w1_Vane.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 10 -folder ~/PARES-PAU/profiles_pk/ -file 'profile_pk_LrM2_w1_Vane.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 10 -folder ~/PARES-PAU/profiles_pk/ -file 'profile_pk_Lrm2_w1_Vane.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 10 -folder ~/PARES-PAU/profiles_pk/ -file 'profile_pk_zm_w1_Vane.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 10 -folder ~/PARES-PAU/profiles_pk/ -file 'profile_pk_zM_w1_Vane.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 10 -folder ~/PARES-PAU/profiles_pk/ -file 'profile_pk_MM_w1_Vane.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 10 -folder ~/PARES-PAU/profiles_pk/ -file 'profile_pk_Mm_w1_Vane.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 10 -folder ~/PARES-PAU/profiles_pk/ -file 'profile_pk_red_w1_Vane.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 10 -folder ~/PARES-PAU/profiles_pk/ -file 'profile_pk_blue_w1_Vane.fits'  -RIN 300 -ROUT 10000 &
wait
