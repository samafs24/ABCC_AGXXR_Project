#!/bin/bash
# Description: Script to identify the links to the FASTQ files from the experiment SDRF file,
#              then download them to the appropriate sample directory (these directories are already created and organized).

echo "Starting FASTQ download process..."

# Set the base directory
baseDir="/scratch/grp/msc_appbio/ABCC_AGXX_Project"
# Set the experiment SDRF file
SDRF="$baseDir/data/experiment_metadata/E-MTAB-7074.sdrf.txt"

echo "Extracting FASTQ file links from the experiment SDRF file: $SDRF"
grep -o 'ftp://[^ ]*\.fastq\.gz' "$SDRF" > "$baseDir/data/experiment_metadata/fastq_links.txt"
echo "FASTQ links successfully extracted."

echo "Beginning download of each FASTQ file to its corresponding directory..."

while IFS= read -r url; do
  echo "Processing URL: $url"
  
  # Extract the sample name from the URL (head to include first match, and tr to exclude newline character)
  sample=$(echo "$url" | grep -o -m1 'ERR[0-9]\{7\}'| head -n 1 | tr -d '\n')
  echo "Sample name extracted: $sample"
  
  # Set the target directory for this sample
 targetDir="$baseDir/data/raw_data/$sample"
    
  # Download the FASTQ file to the target directory
  echo "Downloading file to $targetDir..."
  wget -N -P  "$targetDir" "$url"

done < "$baseDir/data/experiment_metadata/fastq_links.txt"

echo "FASTQ download process complete."



