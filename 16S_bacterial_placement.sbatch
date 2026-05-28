#!/bin/bash
#SBATCH --job-name=silva-bacterial-placement
#SBATCH --partition=msismall
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=falb0011@umn.edu
#SBATCH --output=silva-bacterial-placement-%j.out
#SBATCH --error=silva-bacterial-placement-%j.err

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
cd /users/9/falb0011/FAB2/Summer2022/phylogeneticPlacement/16S_V4_bacteria

# Change to your github script directory (One-time only)
SCRIPT_DIR="/users/9/falb0011/packages/phylogeneticPlacement"

# Change to your output directory (ensure directory exists)
OUTPUT_DIR="/users/9/falb0011/FAB2/Summer2022/phylogeneticPlacement/16S_V4_bacteria"

# Change to your query sequence file
QUERY_FASTA="${OUTPUT_DIR}/FAB2_16S_sequencesQCd.fasta"

# Change to your preferred file name
OUTPUT_NAME="FAB2_2022_16S_phyloPlacement"

SCRATCH_DIR="/scratch.global/falb0011/FAB2_2022_phyloPlacement"


#############################

# Create output folder
mkdir -p "${OUTPUT_DIR}/results"
RESULTS_DIR="${OUTPUT_DIR}/results"

# Create output scratch folder
mkdir -p "$SCRATCH_DIR"

### Load python requirements ###
module load python3
module load infernal/1.1.4
module load cmake
module load conda
source /common/software/install/migrated/anaconda/python3-2020.07-mamba/etc/profile.d/conda.sh
conda init bash
conda activate bacteriaPhyloPlacement

## Specify input tree files ##
REF_TREE_DIR="${SCRIPT_DIR}/refTrees/silva_138.2_nr99_bacteria" # directory to input tree

TREE_FILE="${REF_TREE_DIR}/silva_138.2_nr99_fasttree.nwk"
TREE_LOG="${REF_TREE_DIR}/silva_138.2_nr99_fasttree.log"
ALN_FILE="${REF_TREE_DIR}/silva_138.2_nr99_centroids_masked.fasta"
COLUMN_MASK="${REF_TREE_DIR}/silva_138.2_nr99_centroids_column_mask.npy"

## Specify names of output files ##
QUERY_ALIGNED_RAW="${RESULTS_DIR}/${OUTPUT_NAME}_aligned_raw.pfam"
QUERY_ALIGNED_RAW_FASTA="${RESULTS_DIR}/${OUTPUT_NAME}_aligned_raw.fasta"
QUERY_ALIGNED="${RESULTS_DIR}/${OUTPUT_NAME}_aligned_masked.fasta"


# Step 1: align queries to RF00177.cm to get the same 1533 match columns
echo "Aligning queries to RF00177.cm..."
cmalign --cpu 12 \
    --matchonly \
    --dnaout \
    --outformat Pfam \
    "${REF_TREE_DIR}/RF00177.cm" \
    "$QUERY_FASTA" \
    > "$QUERY_ALIGNED_RAW"

# Convert pfam to fasta afa file
# NOTE: create virtual environment (see README) before running for the first time
echo "Converting pfam alignment to afa (fasta) alignment..."

esl-reformat afa "$QUERY_ALIGNED_RAW" > "$QUERY_ALIGNED_RAW_FASTA"


# Step 2: apply the same column mask used on the reference alignment.
# Reference and query MUST have identical columns for pplacer.
echo "Applying reference column mask..."
python3 "${SCRIPT_DIR}/scripts/mask_alignment.py" \
    "$QUERY_ALIGNED_RAW_FASTA" \
    "$QUERY_ALIGNED" \
    --apply-mask "$COLUMN_MASK"

# Step 2b: remove sequences that are all gaps after masking.
# These arise when cmalign places bases entirely into insert states (--matchonly
# discards them), leaving nothing after the column mask is applied. pplacer
# cannot place them and will crash with "doesn't overlap any reference sequences".
QUERY_ALIGNED_FILTERED="${RESULTS_DIR}/${OUTPUT_NAME}_aligned_masked_filtered.fasta"
export QUERY_ALIGNED QUERY_ALIGNED_FILTERED
echo "Filtering all-gap sequences..."
python3 - <<'EOF'
import os, sys
input_path  = os.environ["QUERY_ALIGNED"]
output_path = os.environ["QUERY_ALIGNED_FILTERED"]
kept = skipped = 0
with open(input_path) as fin, open(output_path, "w") as fout:
    name, seq = "", ""
    for line in fin:
        if line.startswith(">"):
            if name:
                if seq.replace("-", ""):
                    fout.write(f"{name}\n{seq}\n")
                    kept += 1
                else:
                    print(f"  Removed (all gaps): {name.strip()}", file=sys.stderr)
                    skipped += 1
            name, seq = line, ""
        else:
            seq += line.strip()
    if name:
        if seq.replace("-", ""):
            fout.write(f"{name}\n{seq}\n")
            kept += 1
        else:
            print(f"  Removed (all gaps): {name.strip()}", file=sys.stderr)
            skipped += 1
print(f"Kept {kept} sequences, removed {skipped} all-gap sequences.")
EOF
export QUERY_ALIGNED="$QUERY_ALIGNED_FILTERED"

# Step 3: place with pplacer-BSCAMPP
echo "Running pplacer-BSCAMPP..."
t=32

scampp \
    --placement-method pplacer \
    -i "$TREE_LOG" \
    -t "$TREE_FILE" \
    -a "$ALN_FILE" \
    -q "$QUERY_ALIGNED" \
    -d "$SCRATCH_DIR" \
    -o "$OUTPUT_NAME" \
    --threads 16

# Copy the final jplace back to the results directory:
cp "${SCRATCH_DIR}/${OUTPUT_NAME}.jplace" "$RESULTS_DIR/"

echo ""
echo "Placement complete: ${RESULTS_DIR}/${OUTPUT_NAME}.jplace"

# Step 4: construct a tree based on the phylogenetic placements
echo "Running gappa to construct a tree"

gappa examine graft \
    --jplace-path "${RESULTS_DIR}/${OUTPUT_NAME}.jplace" \
    --fully-resolve \
    --out-dir results/

echo "Finished!"