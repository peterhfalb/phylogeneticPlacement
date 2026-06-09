# Phylogenetic placement of amplicon sequences
more to come...

## For Placement of 16S Bacterial amplicon sequences:

First, create a hmmer virtual environment:

```bash
 mamba create -n bacteriaPhyloPlacement -c bioconda -c conda-forge epa-ng hmmer infernal scikit-bio numpy cmake pplacer gappa

conda activate bacteriaPhyloPlacement
pip install bscampp

```

Second, update the 16S_bacterial_placement.sbatch with your file paths, then run it:

```bash
sbatch 16S_bacterial_placement.sbatch
```

Finally, run the Rscript 16S_bacterial_pruneTree.R locally, to prune the tree to only include your query sequences.


## For placement of 18S Protist amplicon sequences:

First create a virtual environment:

```bash

conda create -y -n eukaryotePhyloPlacement -c bioconda -c conda-forge numpy epa-ng infernal gappa

```