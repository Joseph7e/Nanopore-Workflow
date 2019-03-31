#!/usr/bin/sh
#SBATCH --ntasks=1
#SBATCH --mem=409600
#SBATCH --job-name="canu"
#SBATCH --output="canu"

module purge
module load linuxbrew/colsa


reads=$1
sample=$2
genome_size=1000m

canu \
 -p $sample -d $sample.canu \
 genomeSize=$genome_size \
 -nanopore-raw $reads \
 gnuplotTested=true \
 # error with gnuplot so this will skip that step
 useGrid=false\
 gnuplot=/usr/bin/gnuplot\
 ovsMethod=sequential # required with useGrid=false
