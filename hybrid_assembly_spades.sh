#!/usr/bin/sh
#SBATCH --ntasks=1
#SBATCH --mem=409600
#SBATCH --exclude=node117,node118
#SBATCH --job-name="hybrid-fly-only"
#SBATCH --output="hybrid-fly-only"

module purge
module load linuxbrew/colsa

forward=$1
reverse=$2
unpaired=$3
nanopore=$4
sample=$5


spades.py -t 24 -1 $forward -2 $reverse -s $unpaired --nanopore $nanopore  -o spades-hybrid-only-$sample -m 400 --only-assembler
