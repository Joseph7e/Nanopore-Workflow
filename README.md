# Nanopore Workflows
Various commands for handling Nanopore data.
![alt text](sequencing-animated.gif)


# Table of Contents
* [Overview](https://github.com/Joseph7e/Nanopore-Workflow#Overview)  
    * [Basecalling](https://github.com/Joseph7e/Nanopore-Workflow#Basecalling)
    * [Read Processing](https://github.com/Joseph7e/Nanopore-Workflow#Assessing-and-filtering-Nanopore-Data)  
       * [Assessment of data](https://github.com/Joseph7e/Nanopore-Workflow#Accessing-and-filtering-Nanopore-Data)  
       * [Adapter Trimming](https://github.com/Joseph7e/Nanopore-Workflow#Trim-adapters-with-porechop) 
       * [Filter Reads](https://github.com/Joseph7e/Nanopore-Workflow#Accessing-and-filtering-Nanopore-Data)  
    * [Nanopore-only Assembly](https://github.com/Joseph7e/Nanopore-Workflow#Assessing-and-filtering-Nanopore-Data)  
       * [Canu](https://github.com/Joseph7e/Nanopore-Workflow#Accessing-and-filtering-Nanopore-Data)  
       * [Miniasm](https://github.com/Joseph7e/Nanopore-Workflow#Trim-adapters-with-porechop) 
    * [Assembly Polishing](https://github.com/Joseph7e/Nanopore-Workflow#Assessing-and-filtering-Nanopore-Data)

# Overview

How data was produced etc.  
We will use Illumina data in some instances throughout this tutorial (hybrid-assembly, polishing, and assessment). You'll want to run adapter trimming on the illumina reads prior to using them in any instance. IN addition, you will need an illumina-only assembly for the quality assessment of the nanopore assemblies. FOllow my main genome assembly tutorial to produce a high quality assembly and to dtermine information such as insert size and estimated genome size.  

## Example data
 Here we provide some typical test data for nanopore analysis, lambda. Other examples datasets can be found in the SRA, see my other turorials to download this type of data.
  https://www.ncbi.nlm.nih.gov/bioproject/PRJNA477342
  
 
 ```
 # download a lambda dataset
 wget https://www.dropbox.com/s/eml7z2d82n3k8lq/lambda.tar.gz?dl=0 -O lambda.tar.gz
```


## Basecalling

Basecalling is typically run automatically on the sequencing instrument. FOr example the GirdIOn will run the Guppy basecaller as soon as the fast5s are produced. Note that the Guppy installation prodices corrrupted fastqs right now, fix that. The minion may use a different one. 
  
  
Manual: https://denbi-nanopore-training-course.readthedocs.io/en/latest/basecalling/basecalling.html
Alternative Tools: Albacore, DeepNano-blitz, minKNOW, Chiron, Bonito
Comparison of base callers: https://github.com/rrwick/Basecalling-comparison

```
# basecalling with guppy
guppy_basecaller -i <inputdir> -s <output_dir> --flowcell FLO-MIN106 --kit SQK-LSK109 –fast5_out -r -t 15
```


### Repair corrupted read files produced with guppy


The fastq files need to be decompressed for both methods. For my custom workflow, the column for the sort program must match up with the order the fastq was produced (the numbers in the filename). Be sure to test that this works properly.

#### Joes' custom method.
```
cd <fastq_directory>
ls *.fastq | sort -t'_' -k2 -n | xargs cat - > ../raw_reads.fastq
/mnt/lustre/hcgs/joseph7e/scripts/nanopore_fix_fastq.py <fastq_from_above> > <fixed.fastq>
```

#### Nanopore community way
https://community.nanoporetech.com/posts/fastq-errors-on-gridion-an
```
pip install ont-fastq-deconcatenate
apt-get update && apt-get install python3-pip
pip3 install ont-fastq-deconcatenate
fix_concatenated_fastqs -i <path_to_folder_of_fastqs>
```

### Trim adapters with porechop
The porechop "check_reads" option removes the need to specify adapters. It will automatically check and detemrine which ones to remove.
```
porechop --check_reads 1000 -i raw_reads.fastq -o adapter_trimmed.fastq
```

## Filter reads with filtlong
This step is optional. I did not run it through my first attempts. 
```
filtlong --min_mean_q 80 --min_length 2000 <adapter_trimmed.fastq> > filtered.fq
```

## Assessment with Nanoplot
https://github.com/wdecoster/NanoPlot
I usually run this on the raw reads and after any adapter/quality trimming. Run time ~ 2 hrs per 10 GB
```
NanoPlot --fastq <nanopore.fastq> --threads 24 -o <output-dir>
```


# Nanopore only assembly

## Canu

Reads > 1kb an genome size of 130MB

De Bruijn graph contigs were generated with Platanus

```
canu -d canu-assembly -p filt genomeSize=3.5m -nanopore-raw nanopore-reads.fastq

```

## Miniasm

```
mkdir miniasm_assembly
cd miniasm_assembly
sbatch ~/nanopore_assemble.sh ../adapter_trimmed_reads.fastq
```

# Illumina and Nanopore Hybrid Assembly

### Hybrid Assembly w/ spades
https://www.ncbi.nlm.nih.gov/pubmed/26589280
Run this as a normal spades job but specify the --nanopore reads. Note that with low coverage data I found that the nanopore data does not significantly improve the assembly.
```
sbatch /mnt/lustre/hcgs/joseph7e/scripts/GENOME_ASSEMBLY/hybrid_assembly_spades.sh <illumna-forward> <illumina-reverse> <illumina-unpaired> <nanopore-reads> <sample-name>
```
### Hybrid Assembly w/ Masurca


# Genome Assembly Polishing

## Pilon

## Racon

Manuscript:  
Tutoral: https://denbi-nanopore-training-course.readthedocs.io/en/latest/polishing/medaka/racon.html  

Racon polishes a nanopore-only assembly with nanopore read, or illumina data.  
The medaka documentation advises to do four rounds with racon before polishing with medaka since medaka has been trained with racon polished assemblies. We are only doing one round here.

```
# index genome
bwa index <genome.fasta>

# map ont reads to assembly
bwa mem -t 24 -x ont2d <genome.fasta> <nanopore_reads.fastq> > mapping-filteredONT.sam

# polish with racon and produce new consensus sequence
racon -m 8 -x -6 -g -8 -w 500 -t 24 <nanopore_reads.fastq> mapping-filteredONT.sam <genome.fasta> > racon.fasta
```


## Scaffolding w/ LINKS
Note that this process requires a ton of memory. Maybe use a high memory node or used a reduced set of reads.
```
LINKS -f hybrid_assembly_fixed.fasta  -s <txt_file_with_nanopore_read_paths> -b <output_base>
```
```
sbatch ~/nanopore_LINKS.sh <assembly> <nanopore_reads>
```
## Assembly Assessment Scripts
Quast for contiguity, BUSCO for completness, BWA for quality and correctness.
```
quast.py miniasm.fasta
~/scripts/quality_check__genome_busco.sh <assembly.fasta>
```


## Assembly/Read polishing w/ pilon
#### Step 1: Map Illumina reads to Assembly/Read FASTA.
The script below outputs a lot of extra (useful) data. We only need the mapping file.
```
sbatch ~/scripts/bwa_index_and_mapV2.sh <reference.fasta> <forward_reads> <reverse_reads> <sample_name>
mkdir assembly_ploshing && cd assembly_polishing
mv ../bwa_mapping*/sorted_mapped.bam* ./
```

#### Step 2: Prepare genome chunks for array job.
```
grep ">" <assembly.fasta> \
| tr -d ">" \
| shuf \
| split -d -l 200 - genomechunk.
rename genomechunk.0 genomechunk. genomechunk.0*
mkdir pilon
mv genomechunk* pilon/
```

#### Step 3: Edit pilon script and run pilon
At this point you should be in a directory with four files. One mapping file with its index, a directory named pilon (which contains all genome chunk data), and a copy or symlink of your reference assembly.
```
# copy over generic slurm script
cp ~/scripts/nanopore_pilon.slurm pilon.slurm
```
You need to edit a few things in this script.
#SBATCH --array=0-X%8 (change X to equal the number of genome chunks found in pilon dir.
genome="YOUR GENOME FILE HERE"
--frags "NAME OF MAPPING FILE HERE"

```
sbatch ./pilon.slurm
```

#### Step 4: View logs and concatenate polished genome
```
mkdir logs_pilon
mv *.log logs_pilon
cat pilon/*.fasta > polished_genome.fasta
```

#### Step 5: Rinse and repeat
Repeat the entire process on the newly polished genome. Then again and agin, until you're happy.




step 1.) map the reads to the assembly

