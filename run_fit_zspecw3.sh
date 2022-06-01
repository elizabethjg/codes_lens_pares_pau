#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

python -u fit_pairs_profile_withcov.py -ncores 4 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_w3__zspec_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile_withcov.py -ncores 4 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_LrM_w3__zspec_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile_withcov.py -ncores 4 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_Lrm_w3__zspec_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile_withcov.py -ncores 4 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_zM_w3__zspec_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile_withcov.py -ncores 4 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_zm_w3__zspec_best.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile_withcov.py -ncores 4 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_w3__zspec.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile_withcov.py -ncores 4 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_LrM_w3__zspec.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile_withcov.py -ncores 4 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_Lrm_w3__zspec.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile_withcov.py -ncores 4 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_zM_w3__zspec.fits'  -RIN 300 -ROUT 10000 &
python -u fit_pairs_profile_withcov.py -ncores 4 -folder ~/PARES-PAU/profiles_new/ -file 'profile_wc_zm_w3__zspec.fits'  -RIN 300 -ROUT 10000 &
wait
