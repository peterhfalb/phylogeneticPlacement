#!/bin/bash
#SBATCH -job-name=silva-bacterial-placement
#SBATCH -partition=msismall
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=falb0011@umn.edu
#SBATCH --output=silva-bacterial-placement-%.out
#SBATCH --error=silva-bacterial-placement-%.err

# Place query 16S sequences into the SILVA 138.2 NR99 bacterial reference tree
# using epa-ng - SCAMPP.
#
# Usage:
# 1. Create a new folder on the cluster for your project, which includes a fasta file with your query sequences
# 2. Copy the SBATCH file over to the folder for your project
# 3. Update the file paths below to match your setup
# 4. Run the script with sbatch path/to/16S_bacterial_placement.sbatch

set -e


### FILE PATHS TO CHANGE ###
##         ↓↓↓↓↓          ##

# Change to your working directory
cd /users/9/falb0011/YOURWORKINGDIRECTORY

# Change to your github script directory (One-time only)
SCRIPT_DIR="/users/9/falb0011/packages/phylogeneticPlacement"

# Change to your query sequence file
QUERY_FASTA="${SCRIPT_DIR}/YOUR_QUERY_SEQUENCES"

# Change to your output directory (ensure directory exists)
OUTPUT_DIR="/path/to/output/directory"

# Change to your preferred file name
OUTPUT_NAME="nameOfOutputFiles"

#############################

# Create output folder
if [ ! -d "${OUTPUT_DIR}/results"]; then
    mkdir -p "${OUTPUT_DIR}/results"
fi

RESULTS_DIR="${OUTPUT_DIR}/results"

### Load python requirements ###
module load python3

if [ ! -d "${OUTPUT_DIR}/venv"]; then # check if virtual environment exists, install if not
    python3 -m venv "${OUTPUT_DIR}/venv"
    source "${OUTPUT_DIR}/venv/bin/activate"
    pip install -U pip
    pip install -r scikit-bio>=0.5.7 numpy>=1.19.0 bscampp
else
    source "${OUTPUT_DIR}/venv/bin/activate"
fi

## Specify input tree files ##
REF_TREE_DIR="${SCRIPT_DIR}/phyloTrees/silva_138.2_nr99_bacteria_tree" # directory to input tree

TREE_FILE="${REF_TREE_DIR}/silva_138.2_nr99_fasttree.nwk"
TREE_LOG="${REF_TREE_DIR}/silva_138.2_nr99_fasttree.log"
ALN_FILE="${REF_TREE_DIR}/silva_138.2_nr99_centroids_masked.fasta"
COLUMN_MASK="${REF_TREE_DIR}/silva_138.2_nr99_centroids_column_mask.npy"

## Specify names of output files ##
QUERY_ALIGNED_RAW="${OUTPUT_DIR}/${OUTPUT_NAME}_aligned_raw.fasta"
QUERY_ALIGNED="${OUTPUT_DIR}/${OUTPUT_NAME}_aligned_masked.fasta"

module load infernal/1.1.4
module load cmake

# Step 1: align queries to RF00177.cm to get the same 1533 match columns
echo "Aligning queries to RF00177.cm..."
cmalign --cpu 12 \
    --matchonly \
    --dnaout \
    --outformat AFA \
    "${REF_TREE_DIR}/RF00177.cm" \
    "$QUERY_FASTA" \
    > "$QUERY_ALIGNED_RAW"

# Step 2: apply the same column mask used on the reference alignment.
# Reference and query MUST have identical columns for pplacer.
echo "Applying reference column mask..."
python3 "${SCRIPT_DIR}/scripts/mask_alignment.py" \
    "$QUERY_ALIGNED_RAW" \
    "$QUERY_ALIGNED" \
    --apply-mask "$COLUMN_MASK"

# Step 3: place with pplacer-BSCAMPP
echo "Running pplacer-BSCAMPP..."
t=32

scampp \
    --placement-method epa-ng \
    -i "$TREE_LOG" \
    -t "$TREE_FILE" \
    -a "$ALN_FILE" \
    -q "$QUERY_ALIGNED" \
    -d "$OUTPUT_DIR" \
    -o "$OUTPUT_NAME" \

echo ""
echo "Placement complete: ${OUTPUT_DIR}/${OUTPUT_NAME}.jplace"