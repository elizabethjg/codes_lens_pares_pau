#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

##python -u forpair_profile.py -pcat 'Total' -sample 'pk'      -nbins 15 -RIN 100 -ROUT 10000 -ncores 128  -fields w1w3
##python -u forpair_profile.py -pcat 'Total' -sample 'pk_LrM2' -nbins 15 -RIN 100 -ROUT 10000 -ncores 128 -Lratio_min 0.5  -fields w1w3
##python -u forpair_profile.py -pcat 'Total' -sample 'pk_Lrm2' -nbins 15 -RIN 100 -ROUT 10000 -ncores 128 -Lratio_max 0.5  -fields w1w3
##python -u forpair_profile.py -pcat 'Total' -sample 'pk_zM'   -nbins 15 -RIN 100 -ROUT 10000 -ncores 128 -z_min 0.4       -fields w1w3
##python -u forpair_profile.py -pcat 'Total' -sample 'pk_zm'   -nbins 15 -RIN 100 -ROUT 10000 -ncores 128 -z_max 0.4       -fields w1w3
##python -u forpair_profile.py -pcat 'Total' -sample 'pk_red'  -nbins 15 -RIN 100 -ROUT 10000 -ncores 128 -colorcut 'red'  -fields w1w3
##python -u forpair_profile.py -pcat 'Total' -sample 'pk_blue' -nbins 15 -RIN 100 -ROUT 10000 -ncores 128 -colorcut 'blue' -fields w1w3
##python -u forpair_profile.py -pcat 'Total' -sample 'pk_MM2'   -nbins 15 -RIN 100 -ROUT 10000 -ncores 128 -Mtot_min -22.5  -fields w1w3
##python -u forpair_profile.py -pcat 'Total' -sample 'pk_Mm2'   -nbins 15 -RIN 100 -ROUT 10000 -ncores 128 -Mtot_max -22.5  -fields w1w3
##python -u forpair_profile.py -pcat 'Gold' -sample  'pk'      -nbins 15 -RIN 100 -ROUT 10000 -ncores 128  -fields w1w3
##python -u forpair_profile.py -pcat 'Gold' -sample  'pk_LrM2' -nbins 15 -RIN 100 -ROUT 10000 -ncores 128 -Lratio_min 0.5  -fields w1w3
##python -u forpair_profile.py -pcat 'Gold' -sample  'pk_Lrm2' -nbins 15 -RIN 100 -ROUT 10000 -ncores 128 -Lratio_max 0.5  -fields w1w3
##python -u forpair_profile.py -pcat 'Gold' -sample  'pk_zM'   -nbins 15 -RIN 100 -ROUT 10000 -ncores 128 -z_min 0.4       -fields w1w3
##python -u forpair_profile.py -pcat 'Gold' -sample  'pk_zm'   -nbins 15 -RIN 100 -ROUT 10000 -ncores 128 -z_max 0.4       -fields w1w3
##python -u forpair_profile.py -pcat 'Gold' -sample  'pk_red'  -nbins 15 -RIN 100 -ROUT 10000 -ncores 128 -colorcut 'red'  -fields w1w3
##python -u forpair_profile.py -pcat 'Gold' -sample  'pk_blue' -nbins 15 -RIN 100 -ROUT 10000 -ncores 128 -colorcut 'blue' -fields w1w3
##python -u forpair_profile.py -pcat 'Gold' -sample  'pk_MM2'   -nbins 15 -RIN 100 -ROUT 10000 -ncores 128 -Mtot_min -22.5  -fields w1w3
##python -u forpair_profile.py -pcat 'Gold' -sample  'pk_Mm2'   -nbins 15 -RIN 100 -ROUT 10000 -ncores 128 -Mtot_max -22.5  -fields w1w3
python -u forpair_profile.py -pcat 'deepz' -sample  'pk'      -nbins 15 -RIN 100 -ROUT 10000 -ncores 30  -fields w1w3
python -u forpair_profile.py -pcat 'deepz' -sample  'pk_LrM2' -nbins 15 -RIN 100 -ROUT 10000 -ncores 30 -Lratio_min 0.5  -fields w1w3
python -u forpair_profile.py -pcat 'deepz' -sample  'pk_Lrm2' -nbins 15 -RIN 100 -ROUT 10000 -ncores 30 -Lratio_max 0.5  -fields w1w3
python -u forpair_profile.py -pcat 'deepz' -sample  'pk_zM'   -nbins 15 -RIN 100 -ROUT 10000 -ncores 30 -z_min 0.4       -fields w1w3
python -u forpair_profile.py -pcat 'deepz' -sample  'pk_zm'   -nbins 15 -RIN 100 -ROUT 10000 -ncores 30 -z_max 0.4       -fields w1w3
python -u forpair_profile.py -pcat 'deepz' -sample  'pk_red'  -nbins 15 -RIN 100 -ROUT 10000 -ncores 30 -colorcut 'red'  -fields w1w3
python -u forpair_profile.py -pcat 'deepz' -sample  'pk_blue' -nbins 15 -RIN 100 -ROUT 10000 -ncores 30 -colorcut 'blue' -fields w1w3
python -u forpair_profile.py -pcat 'deepz' -sample  'pk_MM2'   -nbins 15 -RIN 100 -ROUT 10000 -ncores 30 -Mtot_min -22.5  -fields w1w3
python -u forpair_profile.py -pcat 'deepz' -sample  'pk_Mm2'   -nbins 15 -RIN 100 -ROUT 10000 -ncores 30 -Mtot_max -22.5  -fields w1w3




