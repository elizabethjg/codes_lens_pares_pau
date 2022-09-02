#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

python -u fit_pairs_profile.py -ncores 7 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_w1w3__vane.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 7 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_LrM2_w1w3__vane.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 7 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_Lrm2_w1w3__vane.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 7 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_zm_w1w3__vane.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 7 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_zM_w1w3__vane.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 7 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_MM_w1w3__vane.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 7 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_Mm_w1w3__vane.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 7 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_red_w1w3__vane.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 7 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_blue_w1w3__vane.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 4 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_M1_w1w3__vane.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 4 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_M2_w1w3__vane.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 4 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_M3_w1w3__vane.fits'  -RIN 300 -ROUT 10000 &
##python -u fit_pairs_profile.py -ncores 4 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_M4_w1w3__vane.fits'  -RIN 300 -ROUT 10000 &
wait
