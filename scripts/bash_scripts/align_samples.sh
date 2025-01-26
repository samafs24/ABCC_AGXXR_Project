#!/bin/bash

# Define paths
baseDir="/scratch/grp/msc_appbio/ABCC_AGXX_Project/data"

referencePrefix="${baseDir}/reference_genome/bowtie2_index/CP000730"

rawDataDir="${baseDir}/raw_data"

processedDataDir="${baseDir}/processed_data"

# Loop through each sample directory
for sampleDir in "$rawDataDir"/*; do
    if [ -d "$sampleDir" ]; then
        # Extract sample name
        sampleName=$(basename "$sampleDir")

        # Define paths for input files and output directory
        R1="${sampleDir}/${sampleName}_1.fastq.gz"
        R2="${sampleDir}/${sampleName}_2.fastq.gz"
        outputDir="${processedDataDir}/${sampleName}"
        outputSAM="${outputDir}/${sampleName}.sam"
	outputBAM="${outputDir}/${sampleName}.bam"

        # Create output directory for this sample
        mkdir -p "$outputDir"

        # Run bowtie2 alignment
        bowtie2 -x "$referencePrefix" -1 "$R1" -2 "$R2" -S "$outputSAM" -p 4
	samtools view -S -b "$outputSAM" > "$outputBAM" # input sam -S output bam -b
    fi
done









