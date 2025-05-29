#!/bin/bash -l
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G
#SBATCH --time=2-00:15:00     # 2 days and 15 minutes
##SBATCH --mail-user=scolando@andrew.cmu.edu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="chipseq"
#SBATCH -p amd # This is the default partition
#SBATCH --output=log_filename.log #if you do not use this, it will create a slurm-[jobid].out as the log

#your other command down here
module load r/4.2.2
module load anaconda3
#bedtools is installed as a module because it is a regular program
module load bedtools

conda activate chipseq
export PATH=$PWD:$PATH
R CMD BATCH HPCrun.R


# #here is how we created the chipseq module
# module load anaconda3
# conda create -n "chipseq" python=3.9
# conda activate chipseq
# conda install pip
# pip install macs2
# pip install MAnorm2-utils

# 
# #and after i'm done creating it type:
# conda deactivate

cp -r /tmpfs/$SLURM_JOB_USER/$SLURM_JOB_ID /rhome/$SLURM_JOB_USER/ChIP-Seq-norm$SLURM_JOB_ID