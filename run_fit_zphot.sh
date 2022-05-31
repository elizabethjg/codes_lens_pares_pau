#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

python -u fit_pairs_profile_withcov.py -ncores 6 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_all__photo_z_2nd_run_mag_i_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile_withcov.py -ncores 6 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_LrM_all__photo_z_2nd_run_mag_i_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile_withcov.py -ncores 6 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_Lrm_all__photo_z_2nd_run_mag_i_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile_withcov.py -ncores 6 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_zM_all__photo_z_2nd_run_mag_i_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile_withcov.py -ncores 6 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_zm_all__photo_z_2nd_run_mag_i_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile_withcov.py -ncores 6 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile_withcov.py -ncores 6 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_LrM_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile_withcov.py -ncores 6 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_Lrm_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile_withcov.py -ncores 6 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_zM_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile_withcov.py -ncores 6 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_zm_all__photo_z_2nd_run_mag_i.fits'  -RIN 300 -ROUT 10000 &
wait
