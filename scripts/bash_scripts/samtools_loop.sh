#!/bin/bash

# Define the sample list
samples=('ERR2713020' 'ERR2713021' 'ERR2713022' 'ERR2713023' 'ERR2713024' 'ERR2713025')

# Loop through each sample
for sample in "${samples[@]}"; do
  # Align the sample using bowtie2
  bowtie2 -x data/reference_genome/CP000730 \
    -1 "data/raw_data/${sample}_1.fastq" \
    -2 "data/raw_data/${sample}_2.fastq" \
    -S "data/processed_data/sam_files/${sample}.sam" \
    -p 4
  
  # Convert SAM to BAM using samtools
  samtools view -S -b "data/processed_data/sam_files/${sample}.sam" > "data/processed_data/bam_files/${sample}.bam"
  
  # Sort the BAM file
  samtools sort "data/processed_data/bam_files/${sample}.bam" -o "data/processed_data/sorted_bam_files/sorted_${sample}.bam"
  
  # Index the sorted BAM file
  samtools index "data/processed_data/sorted_bam_files/sorted_${sample}.bam"
  
  # Run featureCounts
  featureCounts -t CDS -p -g gene_id -F GTF \
    -a data/reference_genome/CP000730.gtf \
    "data/processed_data/sorted_bam_files/sorted_${sample}.bam" \
    -o "data/processed_data/count_files/${sample}_counts.txt"
done

# Prepare a list of BAM files for all samples
bam_files=""
for sample in "${samples[@]}"; do
  bam_files+="data/processed_data/sorted_bam_files/sorted_${sample}.bam "
done

# Run featureCounts for all BAM files and generate a single output file
featureCounts -t CDS -p -g gene_id -F GTF \
  -a data/reference_genome/CP000730.gtf \
  $bam_files \
  -o data/processed_data/count_files/all_samples_counts.txt
