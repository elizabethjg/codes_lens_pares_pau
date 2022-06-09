#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv
python -u fit_pairs_profile.py -ncores 6 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_all__zspec_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 6 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_LrM_all__zspec_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 6 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_Lrm_all__zspec_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 6 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_zM_all__zspec_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 6 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_zm_all__zspec_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 6 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_all__zspec.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 6 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_LrM_all__zspec.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 6 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_Lrm_all__zspec.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 6 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_zM_all__zspec.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile.py -ncores 6 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_zm_all__zspec.fits'  -RIN 300 -ROUT 10000 &
wait
