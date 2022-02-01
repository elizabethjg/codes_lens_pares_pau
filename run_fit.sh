#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

python -u fit_pairs_profile.py -ncores 56 -folder ~/PARES-PAU/profiles/ -file 'profile_all_100_10000_photo_z_2nd_run_mag_i_mask.fits'  -RIN 200 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 56 -folder ~/PARES-PAU/profiles/ -file 'profile_all_100_10000_zspec_mask.fits'  -RIN 200 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 56 -folder ~/PARES-PAU/profiles/ -file 'profile_all_200_10000_photo_z_2nd_run_mag_i_mask.fits'  -RIN 200 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 56 -folder ~/PARES-PAU/profiles/ -file 'profile_all_200_10000_zspec_mask.fits' -RIN 200 -ROUT 10000 &
wait
