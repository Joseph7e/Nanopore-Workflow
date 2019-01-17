# Repair corrupted read files prodiced with guppy
### Joes Way
```
cd <fastq_directory>
ls *.fastq | sort -t'_' -k2 -n | xargs cat - > ../raw_reads.fastq
/mnt/lustre/hcgs/joseph7e/scripts/nanopore_fix_fastq.py <fastq_from_above> > <fixed.fastq>
```

### Nanopore community way
https://community.nanoporetech.com/posts/fastq-errors-on-gridion-an
```
pip install ont-fastq-deconcatenate
apt-get update && apt-get install python3-pip
pip3 install ont-fastq-deconcatenate
fix_concatenated_fastqs -i <path_to_folder_of_fastqs>
```

# trim adapters
```
porechop --check_reads 1000 -i raw_reads.fastq -o adapter_trimmed.fastq
```
# Optional: filter reads based on length and quality
filtlong --min_mean_q 80 --min_length 2000 <adapter_trimmed.fastq> > filtered.fq

# Examine reads w/ NanoPlot
```
NanoPlot --fastq <nanopore.fastq> --threads 24 -o <output-dir>
```
# Nanopore only assembly w/ Miniasm + Minimap
```
mkdir miniasm_assembly
cd miniasm_assembly
sbatch ~/nanopore_assemble.sh ../adapter_trimmed_reads.fastq
```

# Hybrid Assembly w/ spades

# Hybrid Assembly w/ Masurca

# Scaffolding w/ LINKS
Note that this process requires a ton of memory. Maybe use a high memory node or used a reduced set of reads.
```
LINKS -f hybrid_assembly_fixed.fasta  -s <txt_file_with_nanopore_read_paths> -b <output_base>
```
```
sbatch ~/nanopore_LINKS.sh <assembly> <nanopore_reads>
```
# assess
```
quast.py miniasm.fasta
~/scripts/quality_check__genome_busco.sh <assembly.fasta>
```
# Assembly/Read polishing w/ pilon




