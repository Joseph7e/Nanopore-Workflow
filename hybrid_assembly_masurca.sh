#!/usr/bin/sh
#SBATCH --ntasks=1
#SBATCH --mem=409600
#SBATCH --job-name="masurca_E02"
#SBATCH --output="masurca_E02"

module purge
module load linuxbrew/colsa

#sample="Sample_184615539"


#export AUGUSTUS_CONFIG_PATH="/mnt/lustre/software/linuxbrew/colsa/Cellar/augustus/3.2.2_2/libexec/config/"
#


#### concatenate files
sample=$1
raw_forward=$2 # USE ABSOLUTE PATHS!
raw_reverse=$3
insert_metrics=$4
genome_cov_file=$5 # concoct_input_table.tsv
optional_filter_list=$6
nanopore=$7


cwd=$(pwd)
echo $cwd

echo Gathering data for config file
# calculate insert sizes
echo grep -A 1 MEDIAN_INSERT $insert_metrics | grep -v 'MEDIAN_INSERT' | awk -F'\t' '{print $6}'
INSERT_SIZE="$(grep -A 1 MEDIAN_INSERT $insert_metrics | grep -v 'MEDIAN_INSERT' | awk -F'\t' '{print $6}')"
INSERT_STDV="$(grep -A 1 MEDIAN_INSERT $insert_metrics | grep -v 'MEDIAN_INSERT' | awk -F'\t' '{print $7}')"

echo INSERT SIZE $INSERT_SIZE
echo INSERT STDV $INSERT_STDV


# calculate genome and coverage estimates
#echo /mnt/lustre/hcgs/joseph7e/scripts/GENOME_ASSEMBLY/length_and_avg_cov_from_concoct_input_table.py $genome_cov_file $optional_filter_list
string="$(/mnt/lustre/hcgs/joseph7e/scripts/GENOME_ASSEMBLY/length_and_avg_cov_from_concoct_input_table.py $genome_cov_file $optional_filter_list)"
GENOME_SIZE=$(echo $string | cut -f1 -d' ')
EST_COVERAGE=$(echo $string | cut -f2 -d' ')

mkdir masurca_assembly_$sample && cd masurca_assembly_$sample
cd masurca_assembly_$sample

# make a copy of generic masurca configs
cp /mnt/lustre/hcgs/joseph7e/scripts/masurca_example_config.txt ./

# add reads into config file
sed -i "s:FORWARD_READS:../$raw_forward:g" masurca_example_config.txt
sed -i "s:REVERSE_READS:../$raw_reverse:g" masurca_example_config.txt

# adjust insert sizes in config file
sed -i "s:INSERT_SIZE_AVERAGE:$INSERT_SIZE:g" masurca_example_config.txt
sed -i "s:INSERT_SIZE_STDV:$INSERT_STDV:g" masurca_example_config.txt


# adjust JF hash value, use genome_size * estimated_coverage
echo genome_size $GENOME_SIZE
echo estimated_coverage $EST_COVERAGE
p=$(echo "$GENOME_SIZE * $EST_COVERAGE" | bc)
# round the number to get final value
JF_SIZE=$(echo ${p%%.*})
echo JellyFish_Hash_Size $JF_SIZE
sed -i "s:GENOME_JELLYFISH_HASH:$JF_SIZE:g" masurca_example_config.txt


echo running MASURCA:

module purge
module load anaconda/colsa
source activate masurca-3.2.6

# run masurca
masurca masurca_example_config.txt
srun ./assemble.sh


#    echo srun spades.py -1 $forward -2 $reverse -s $unpaired -t 24 -o $dir/spades_assembly
#   # srun spades.py -1 $forward -2 $reverse -s $unpaired -t 24 -o $dir/spades_assembly
#    srun spades.py -o $dir/spades_assembly --continue -t 24 -m 400
#cd ../../
#fasta=$dir/masurca_assembly/CA/10-gapclose/genome.ctg.fasta
#scaffold=$dir/masurca_assembly/CA/10-gapclose/genome.scf.fasta

#echo '#'running QUAST:
#echo mkdir QUAST_results
#srun quast.py -o $dir/QUAST_results/quast_masurca_$sample $fasta
#srun quast.py -o $dir/QUAST_results/quast_masurca_scaff_$sample $fscaffold
##    run quast.py -o $dir/QUAST_results/quast_$sample $fasta
#
#
#echo '#'running BUSCO:
#echo srun busco -i $fasta -m genome --cpu 24 -o $sample
##srun busco -i $fasta -m genome --cpu 24 -o $sample

