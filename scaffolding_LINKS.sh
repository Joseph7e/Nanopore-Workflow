#!/usr/bin/sh
#SBATCH --ntasks=1
#SBATCH --mem=409600
#SBATCH --job-name="LINKS"
#SBATCH --output="slurm.LINKS"

module purge
module load linuxbrew/colsa


genome=$1
reads=$2
output=LINKS-scaffolding

readlink -f $reads > nanopore_read_locations.txt

LINKS -f $genome -s nanopore_read_locations.txt -b $output
