#!/bin/bash
#SBATCH --job-name=profpairs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=3-00:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


. /etc/profile
source $HOME/.bashrc

##cd /mnt/projects/lensing/HALO_SHAPE/MICEv1.0/codes_HS_micev1.0/

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

conda activate myenv

srun python -u forpair_profile.py -ncores 56 -sample 'w3' -nbins 20 -RIN 100 -ROUT 10000
##srun python -u forpair_profile.py -ncores 56 -sample 'all' -nbins 20 -RIN 100 -ROUT 10000
##srun python -u forpair_profile.py -ncores 56 -sample 'zL' -nbins 20 -RIN 100 -ROUT 10000 -z_min 0. -z_max 0.4
##srun python -u forpair_profile.py -ncores 56 -sample 'zH' -nbins 20 -RIN 100 -ROUT 10000 -z_min 0.4 
##srun python -u forpair_profile.py -ncores 56 -sample 'LratioL' -nbins 20 -RIN 100 -ROUT 10000 -Lratio_min 0. -Lratio_max 0.5
##srun python -u forpair_profile.py -ncores 56 -sample 'LratioH' -nbins 20 -RIN 100 -ROUT 10000 -Lratio_min 0.5 

