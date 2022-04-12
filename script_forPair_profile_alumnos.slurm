#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

#python -u forpair_profile.py -sample 'new' -nbins 20 -RIN 100 -ROUT 10000 -ncores 32 -fields 'all'
#python -u forpair_profile.py -sample 'new' -nbins 20 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w3'
#python -u forpair_profile.py -sample 'new' -nbins 20 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1w2'
python -u forpair_profile.py -sample 'new' -nbins 20 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w2w3'
python -u forpair_profile.py -sample 'new' -nbins 20 -RIN 100 -ROUT 10000 -ncores 32 -fields 'w1w3'
##srun python -u forpair_profile.py -ncores 56 -sample 'all' -nbins 20 -RIN 100 -ROUT 10000
##srun python -u forpair_profile.py -ncores 56 -sample 'zL' -nbins 20 -RIN 100 -ROUT 10000 -z_min 0. -z_max 0.4
##srun python -u forpair_profile.py -ncores 56 -sample 'zH' -nbins 20 -RIN 100 -ROUT 10000 -z_min 0.4 
##srun python -u forpair_profile.py -ncores 56 -sample 'LratioL' -nbins 20 -RIN 100 -ROUT 10000 -Lratio_min 0. -Lratio_max 0.5
##srun python -u forpair_profile.py -ncores 56 -sample 'LratioH' -nbins 20 -RIN 100 -ROUT 10000 -Lratio_min 0.5 

