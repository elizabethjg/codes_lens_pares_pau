#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

python -u fit_pairs_profile.py -ncores 15 -folder ~/PARES-PAU/profiles/ -file 'profile_newall_100_10000_zback_photo_z_2nd_run_mag_i_mask.fits'  -RIN 200 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 15 -folder ~/PARES-PAU/profiles/ -file 'profile_newall_100_10000_zback_zspec_mask.fits'  -RIN 200 -ROUT 10000 &
wait
