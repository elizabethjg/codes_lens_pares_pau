#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

python -u fit_pairs_profile.py -ncores 3 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_LrM_all__photo_z_2nd_run_mag_i_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 3 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_Lrm_all__photo_z_2nd_run_mag_i_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 3 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_zm_all__photo_z_2nd_run_mag_i_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 3 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_zM_all__photo_z_2nd_run_mag_i_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 3 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_LrM_w3__photo_z_2nd_run_mag_i_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 3 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_Lrm_w3__photo_z_2nd_run_mag_i_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 3 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_zm_w3__photo_z_2nd_run_mag_i_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 3 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_zM_w3__photo_z_2nd_run_mag_i_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 3 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 3 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_LrM_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 3 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_Lrm_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 3 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_zm_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 3 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_zM_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 3 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_w3__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 3 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_LrM_w3__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 3 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_Lrm_w3__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 3 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_zm_w3__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 3 -folder ~/PARES-PAU/profiles3/ -file 'profile_mh_zM_w3__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
wait
