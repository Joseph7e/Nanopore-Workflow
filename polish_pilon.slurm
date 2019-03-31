#!/bin/bash

#SBATCH --partition=shared
#SBATCH -J pilon
#SBATCH --cpus-per-task=24
#SBATCH --output pilon.%A_%a.log
#SBATCH --array=0-42%8
##SBATCH --dependency=afterok:57955,57956,57957,57958,57959,57960
#SBATCH --mem 110Gb

#source ~/.bash_profile
#PATH=$PATH:$HOME/Oyster_River_Protocol/software/anaconda/install/bin/
#source activate gapcloser

module load linuxbrew/colsa

genome=miniasm.fasta

#cd /mnt/lustre/macmaneslab/macmanes/pero_genome/

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

pilon --genome $genome \
--frags sorted_mapped.bam \
--output pilon/pilonchunk.$SLURM_ARRAY_TASK_ID \
--fix bases,gaps \
--diploid \
--threads 24 \
--flank 5 \
--verbose \
--mingap 1 \
--nostrays \
--targets pilon/genomechunk.$SLURM_ARRAY_TASK_ID
