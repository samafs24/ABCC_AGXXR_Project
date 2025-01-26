#!/bin/bash

echo "Start Bowtie2-build "

# Define the base directory and file path
baseDir="/scratch/grp/msc_appbio/ABCC_AGXX_Project/data/reference_genome"
referenceGenome="${baseDir}/CP000730.1.fasta"
outputPrefix="${baseDir}/CP000730"

# Run Bowtie2-build
bowtie2-build -f --threads 4 "${baseDir}/CP000730.1.fasta" "${baseDir}/CP000730"

echo "End Bowtie2-build "
