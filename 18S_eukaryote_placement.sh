#!/bin/bash
#SBATCH --job-name=eukaryote-placement
#SBATCH --partition=msismall
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=falb0011@umn.edu
#SBATCH --output=eukaryote-placement-%j.out
#SBATCH --error=eukaryote-placement-%j.err

# Place query 18S SSU amplicon sequences into the Romero et al. SSU+LSU eukaryote
# reference tree using epa-ng.
#
# The reference tree was built from concatenated 18S (SSU) + 28S (LSU) sequences.
# Query sequences (18S only) are aligned to the RF01960 covariance model, then padded
# with dashes for the LSU portion so their column count matches the reference alignment.
#
# Usage:
# 1. Create a new folder on the cluster for your project with your query fasta file
# 2. Copy this script to that folder
# 3. Update the file paths below to match your setup
# 4. Run with: sbatch path/to/18S_eukaryote_placement.sh
#
# Dependencies (must be available in your environment):
#   infernal (cmalign, esl-reformat)
#   epa-ng
#   gappa
#   python3 with numpy

set -e


### FILE PATHS TO CHANGE ###
##         ↓↓↓↓↓          ##

# Change to your working directory
cd /users/9/falb0011/FAB2/Summer2022/phylogeneticPlacement/18S_V4_protist

# Change to your github script directory (one-time clone)
SCRIPT_DIR="/users/9/falb0011/packages/phylogeneticPlacement"

# Change to your output directory (ensure directory exists)
OUTPUT_DIR="/users/9/falb0011/FAB2/Summer2022/phylogeneticPlacement/18S_V4_protist"

# Change to your query 18S sequence file (unaligned fasta)
QUERY_FASTA="${OUTPUT_DIR}/FAB2.18S.protist.sequencesQCd.fasta"

# Change to your preferred output file prefix
OUTPUT_NAME="FAB2_2022_18S_phyloPlacement"

#############################

# Create results directory
mkdir -p "${OUTPUT_DIR}/results"
RESULTS_DIR="${OUTPUT_DIR}/results"

### Load conda and activate environment ###
# One-time setup (run once before first use):
#   conda create -y -n eukaryotePhyloPlacement -c bioconda -c conda-forge \
#       epa-ng infernal gappa numpy hmmer
module load python3/3.12.4_anaconda2024.06-1_libmamba
conda init bash
conda activate eukaryotePhyloPlacement

## Reference tree files (do not change) ##
REF_TREE_DIR="${SCRIPT_DIR}/refTrees/eukaryote_ssu.lsu_tree"
TREE_FILE="${REF_TREE_DIR}/3dom-org.18s28s.cmalign.f.fasta.contree"
ALN_FILE="${REF_TREE_DIR}/3dom-org.18s28s.cmalign.f.fasta"
CM_FILE="${REF_TREE_DIR}/RF01960.cm"

## Intermediate and output file names ##
QUERY_ALIGNED_RAW="${RESULTS_DIR}/${OUTPUT_NAME}_aligned_raw.pfam"
QUERY_ALIGNED_SSU="${RESULTS_DIR}/${OUTPUT_NAME}_aligned_ssu.fasta"
QUERY_ALIGNED="${RESULTS_DIR}/${OUTPUT_NAME}_aligned_padded.fasta"


# Step 1: align queries to RF01960.cm to get the 18S SSU match columns
if [ -f "$QUERY_ALIGNED_RAW" ]; then
    echo "Skipping cmalign — pfam file already exists: $QUERY_ALIGNED_RAW"
else
    echo "Aligning queries to RF01960.cm..."
    cmalign --cpu 12 \
        --matchonly \
        --dnaout \
        --outformat Pfam \
        --mxsize 8096 \
        "$CM_FILE" \
        "$QUERY_FASTA" \
        > "$QUERY_ALIGNED_RAW"
fi

# Convert pfam to fasta
if [ -f "$QUERY_ALIGNED_SSU" ]; then
    echo "Skipping esl-reformat — SSU fasta already exists: $QUERY_ALIGNED_SSU"
else
    echo "Converting pfam alignment to fasta..."
    esl-reformat afa "$QUERY_ALIGNED_RAW" > "$QUERY_ALIGNED_SSU"
fi


# Step 2: pad SSU alignment with dashes for the LSU (28S) portion.
# The reference alignment is SSU+LSU concatenated; query sequences only cover SSU.
# Dashes are appended so the column count matches the reference alignment exactly.
export QUERY_ALIGNED_SSU ALN_FILE QUERY_ALIGNED
if [ -f "$QUERY_ALIGNED" ]; then
    echo "Skipping padding — padded fasta already exists: $QUERY_ALIGNED"
else
    echo "Padding SSU alignment with LSU dashes..."
    python3 - <<'EOF'
import os, sys

ssu_fasta = os.environ["QUERY_ALIGNED_SSU"]
ref_fasta = os.environ["ALN_FILE"]
out_fasta = os.environ["QUERY_ALIGNED"]

# Determine reference alignment total column count
ref_len = None
with open(ref_fasta) as f:
    seq = ""
    for line in f:
        if line.startswith(">"):
            if seq:
                ref_len = len(seq)
                break
        else:
            seq += line.strip()
if ref_len is None:
    sys.exit("ERROR: Could not determine reference alignment length")
print(f"Reference alignment length: {ref_len} columns")

# Read SSU alignment
ssu_seqs = {}
with open(ssu_fasta) as f:
    name = None
    for line in f:
        if line.startswith(">"):
            name = line.strip()
            ssu_seqs[name] = ""
        elif name:
            ssu_seqs[name] += line.strip()

if not ssu_seqs:
    sys.exit("ERROR: No sequences found in SSU alignment")

ssu_len = len(next(iter(ssu_seqs.values())))
lsu_pad = ref_len - ssu_len
print(f"SSU alignment length: {ssu_len} columns (from RF01960 match states)")
print(f"LSU pad length:       {lsu_pad} dashes")

if lsu_pad < 0:
    sys.exit(f"ERROR: SSU alignment ({ssu_len}) is longer than reference ({ref_len})")

dash_pad = "-" * lsu_pad
with open(out_fasta, "w") as f:
    for name, seq in ssu_seqs.items():
        f.write(f"{name}\n{seq}{dash_pad}\n")
print(f"Wrote {len(ssu_seqs)} padded sequences to {out_fasta}")
EOF
fi


# Step 2b: remove sequences that are all gaps in the SSU portion.
# These arise when cmalign places all bases into insert states; epa-ng cannot
# place them and will error.
QUERY_ALIGNED_FILTERED="${RESULTS_DIR}/${OUTPUT_NAME}_aligned_padded_filtered.fasta"
export QUERY_ALIGNED QUERY_ALIGNED_FILTERED
if [ -f "$QUERY_ALIGNED_FILTERED" ]; then
    echo "Skipping gap filter — filtered fasta already exists: $QUERY_ALIGNED_FILTERED"
else
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
fi


# Step 3: place sequences into the reference tree with epa-ng.
# Query sequences are already aligned (same column count as reference), so
# epa-ng uses split mode and skips internal hmmer re-alignment.
JPLACE="${RESULTS_DIR}/${OUTPUT_NAME}.jplace"
if [ -f "$JPLACE" ]; then
    echo "Skipping epa-ng — jplace already exists: $JPLACE"
else
    echo "Running epa-ng..."
    epa-ng \
        --tree "$TREE_FILE" \
        --ref-msa "$ALN_FILE" \
        --query "$QUERY_ALIGNED_FILTERED" \
        --redo \
        --out-dir "$RESULTS_DIR" \
        --prefix "$OUTPUT_NAME" \
        --threads 16
    echo ""
    echo "Placement complete: $JPLACE"
fi


# Step 4: graft query placements onto the reference tree to produce a newick
NEWICK="${RESULTS_DIR}/${OUTPUT_NAME}.newick"
if [ -f "$NEWICK" ]; then
    echo "Skipping gappa — newick already exists: $NEWICK"
else
    echo "Running gappa to construct a grafted tree..."
    gappa examine graft \
        --jplace-path "$JPLACE" \
        --fully-resolve \
        --out-dir "${RESULTS_DIR}/"
fi

echo "Finished! Output tree: $NEWICK"
