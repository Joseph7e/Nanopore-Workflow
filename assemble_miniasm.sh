#!/usr/bin/sh
#SBATCH --ntasks=1
##SBATCH --mem=409600
#SBATCH --job-name="miniasm"
#SBATCH --output="miniasm"

module purge
module load linuxbrew/colsa


minimap2 -x ava-ont $1 $1 | gzip -1 > mapped.paf.gz
miniasm -f $1 mapped.paf.gz > miniasm.gfa
awk '/^S/{print ">"$2"\n"$3}' miniasm.gfa > miniasm.fasta



