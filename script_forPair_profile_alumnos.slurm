#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

python -u forpair_profile.py -sample 'new' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w3' 
python -u forpair_profile.py -sample 'new_LrM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w3' -Lratio_min 0.5
python -u forpair_profile.py -sample 'new_Lrm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w3' -Lratio_max 0.5
python -u forpair_profile.py -sample 'new_zM' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w3' -z_min 0.4
python -u forpair_profile.py -sample 'new_zm' -nbins 15 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w3' -z_max 0.4
