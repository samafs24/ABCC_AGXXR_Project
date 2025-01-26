#!/bin/bash
#SBATCH --job-name=abcc_pipeline         # Job name
#SBATCH --partition=msc_appbio           # Partition (queue) name
#SBATCH --nodes=1                        # Number of nodes
#SBATCH --ntasks=1                       # Number of tasks per node
#SBATCH --cpus-per-task=4                # Number of CPU cores per task
#SBATCH --mem=16G                        # Memory per node
#SBATCH --time=12:00:00                  # Time limit (hh:mm:ss)
#SBATCH --output=abcc_pipeline_%j.log    # Standard output and error log

# Load modules if needed
echo "Loading required modules..."
module load anaconda3
module load bowtie2
module load samtools
module load fastqc
module load subread
echo "Modules loaded successfully."

# Check if the user provided a home path and environemnt argument
if [ $# -lt 2 ]; then
    echo "Error: Insufficient arguments provided."
    echo "Usage: $0 <home_path> <conda_env_name>"
    exit 1
fi

# Read the arguments
HOME_PATH=$1
if [ ! -d "$HOME_PATH" ]; then
    echo "Error: The specified home path '$HOME_PATH' does not exist."
    exit 1
fi

CONDA_ENV=$2
echo "Home path set to: $HOME_PATH"
echo "Conda environment name set to: $CONDA_ENV"

# Go to the home path 
echo "Navigating to the home path..."
cd "$HOME_PATH" || exit
echo "Successfully navigated to $HOME_PATH."

# Create the directory structure
echo "Creating directory structure..."
mkdir -p data
mkdir -p data/reference_genome
mkdir -p data/raw_data
mkdir -p data/experiment_metadata
mkdir -p data/processed_data
mkdir -p data/processed_data/fastqc_results
mkdir -p data/processed_data/sam_files
mkdir -p data/processed_data/bam_files
mkdir -p data/processed_data/sorted_bam_files
mkdir -p data/processed_data/count_files
echo "Directory structure created successfully."

# Set up the Conda environment
echo "Checking Conda version..."
export PATH="/opt/anaconda3/bin:$PATH"
conda --version

# Create and activate the Conda environment
if ! conda info --envs | grep -q "$CONDA_ENV"; then
    echo "Creating Conda environment '$CONDA_ENV'..."
    conda create --name "$CONDA_ENV" -y
    echo "Conda environment '$CONDA_ENV' created."
else 
    echo "Conda environment '$CONDA_ENV' already exists."
fi

echo "Activating Conda environment '$CONDA_ENV'..."
source activate "$CONDA_ENV"
conda install -y bowtie2 samtools fastqc subread
echo "Conda environment '$CONDA_ENV' is ready."

# Download reference genome and annotation files
echo "Downloading reference genome and annotation files..."
cd data/reference_genome || { echo "Error: Failed to navigate to data/reference_genome"; exit 1; }
# From RefSeq (GCF)
#wget -N https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/017/085/GCF_000017085.1_ASM1708v1/GCF_000017085.1_ASM1708v1_genomic.fna.gz -O CP000730.fna.gz
#wget -N https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/017/085/GCF_000017085.1_ASM1708v1/GCF_000017085.1_ASM1708v1_genomic.gtf.gz -O CP000730.gtf.gz
# From GenBank (it seems our original study uses this) 
wget -N https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/085/GCA_000017085.1_ASM1708v1/GCA_000017085.1_ASM1708v1_genomic.fna.gz -O CP000730.fna.gz
wget -N https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/085/GCA_000017085.1_ASM1708v1/GCA_000017085.1_ASM1708v1_genomic.gtf.gz -O CP000730.gtf.gz
gunzip -f *.gz
# mv GCF_000017085.1_ASM1708v1_genomic.fna CP000730.fna
# mv GCF_000017085.1_ASM1708v1_genomic.gtf CP000730.gtf
echo "Reference genome and annotation files downloaded and processed."

# Build the reference index
echo "Building Bowtie2 reference genome index..."
bowtie2-build --threads 4 CP000730.fna CP000730
echo "Reference genome index built successfully."

# Download metadata and raw data
echo "Downloading metadata and raw data..."
cd ../..
wget -O data/experiment_metadata/E-MTAB-7074.sdrf.txt https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7074/E-MTAB-7074.sdrf.txt

fastq_links=($(grep -o 'ftp://[^ ]*\.fastq\.gz' data/experiment_metadata/E-MTAB-7074.sdrf.txt))

# Initialize an empty array for sample names
samples=()

# Download FASTQ files and populate the sample list
echo "Downloading FASTQ files and extracting sample names..."
for url in "${fastq_links[@]}"; do
  echo "Processing URL: $url"
  
  # Extract the sample name from the URL
  sample=$(echo "$url" | grep -o -m1 'ERR[0-9]\{7\}'| head -n 1 | tr -d '\n')
  echo "Sample name extracted: $sample"
  
  # Add the sample name to the list
  samples+=("$sample")
  
  # Download the FASTQ file
  wget -N -P data/raw_data "$url"
done
echo "Metadata and raw data downloaded and uncompressed."
echo "Sample names extracted: ${samples[@]}"

# Run FastQC on the raw data
echo "Running FastQC on raw data..."
fastqc -o data/processed_data/fastqc_results -t 4 data/raw_data/*.fastq
echo "FastQC analysis completed."

# Align reads and process files
echo "Aligning reads and processing files..."
for sample in "${samples[@]}"; do
  echo "Processing sample: $sample"
  echo "Starting alignment with Bowtie2 for sample: $sample..."
  bowtie2 -x data/reference_genome/CP000730 \
    -1 "data/raw_data/${sample}_1.fastq" \
    -2 "data/raw_data/${sample}_2.fastq" \
    -S "data/processed_data/sam_files/${sample}.sam" \
    -p 4
  echo "Bowtie2 alignment completed for sample: $sample."

  echo "Converting SAM to BAM for sample: $sample..."
  samtools view -S -b "data/processed_data/sam_files/${sample}.sam" > "data/processed_data/bam_files/${sample}.bam"
  echo "SAM to BAM conversion completed for sample: $sample."

  # Aleksej said we didn't need this; we might want to uncomment it to visualize our genome in IVG
  #echo "Sorting BAM file for sample: $sample..."
  #samtools sort "data/processed_data/bam_files/${sample}.bam" -o "data/processed_data/sorted_bam_files/sorted_${sample}.bam"
  #echo "BAM sorting completed for sample: $sample."

  #echo "Indexing sorted BAM file for sample: $sample..."
  #samtools index "data/processed_data/sorted_bam_files/sorted_${sample}.bam"
  #echo "Indexing completed for sample: $sample."

  echo "Sample $sample processing completed successfully!"
done

# Run featureCounts
echo "Running featureCounts on aligned files..."
#bam_files=$(ls data/processed_data/sorted_bam_files/sorted_*.bam | tr '\n' ' ')
bam_files=$(ls data/processed_data/bam_files/*.bam | tr '\n' ' ')
featureCounts -t CDS -p -g gene_id -F GTF \
  -a data/reference_genome/CP000730.gtf \
  -o data/processed_data/count_files/all_samples_counts.txt $bam_files
echo "featureCounts completed."

Rscript r_script

# Final message
echo "Pipeline completed successfully!"
