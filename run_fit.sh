#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

python -u fit_pairs_profile.py -ncores 13 -folder ~/PARES-PAU/profiles/ -file 'profile_all_photo_z_2nd_run_mag_i_mask.fits'  -RIN 200 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 13 -folder ~/PARES-PAU/profiles/ -file 'profile_zL_photo_z_2nd_run_mag_i_mask.fits'  -RIN 200 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 13 -folder ~/PARES-PAU/profiles/ -file 'profile_zH_photo_z_2nd_run_mag_i_mask.fits'  -RIN 200 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 13 -folder ~/PARES-PAU/profiles/ -file 'profile_LratioL_photo_z_2nd_run_mag_i_mask.fits'  -RIN 200 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 13 -folder ~/PARES-PAU/profiles/ -file 'profile_LratioH_photo_z_2nd_run_mag_i_mask.fits'  -RIN 200 -ROUT 10000 &
wait
