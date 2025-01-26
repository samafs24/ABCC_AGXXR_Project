#!/bin/bash

# Load FastQC module
module load fastqc

# Define base and results directories
baseDirectory="/scratch/grp/msc_appbio/ABCC_AGXX_Project/data/raw_data"
resultsDirectory="/scratch/grp/msc_appbio/ABCC_AGXX_Project/data/processed_data/fastqc_output"

# Create the results directory if it does not exist
mkdir -p "$resultsDirectory"

# Loop through each sample subdirectory
for sampleDir in "$baseDirectory"/ERR*; do
    # Check if the directory exists and contains .fastq.gz files
    if [ -d "$sampleDir" ]; then
        echo "Running FastQC on $sampleDir"

        # Run FastQC on all .fastq.gz files in the sample directory
        fastqc -o "$resultsDirectory" -t 4 "$sampleDir"/*.fastq.gz
    else
        echo "No .fastq.gz files found in $sampleDir."
    fi
done

echo "End of the pipeline"
