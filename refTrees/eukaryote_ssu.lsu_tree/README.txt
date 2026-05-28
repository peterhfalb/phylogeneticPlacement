
## README.txt ##

# This is a backbone phylogeny created by Miguel Romero Gutierrez
# It is part of the analysis referenced in this preprint: https://doi/org/10/1101/2025.03.26.645542

# It is a SSU-LSU concatenated phylogeny built using clustered SILVA, PR2, long-read amplicon datasets,
# and assembled LSU + SSU contigs from metagenome datasets
# It represents probably the best tree we currently have for the eukaryote 18S/28S genes
# It represents a pretty balanced taxon sampling across the tree and should have good representation within Fungi too
# Sequences in the tree are clustered at 85% similarity, so the primary limitation would be that it does not necessarily represent fine-scale diversity

# Here is what Frederik Schulz said about the tree:

Regarding mapping amplicons to the backbone tree, yes, this is an approach that works well, especially when you have relatively closely related sequences in your alignment/tree.

But if your amplicon is too novel (it's difficult to define a definite threshold), adding the full length sequence might change the tree's topology to place it, but if it's mapped it might end up at an internal node and could be mislabeled. To get reliable results, describe what maps with high confidence, but then you can still discuss sequences that had lower confidence mapping to internal nodes. Miguel has done this for his project and should be able to provide more detailed information.

Miguel will be able to share our 18+28S "conservative" tree + alignment with you, it was reasonably well resolved considering that we haven't used phylogenomics and it has a very good and balanced taxon sampling as it takes into account metagenomic OTUs. These OTUs are only 85% sequence similarity OTUs, so they might not be perfect if you need species-level classification. Miguel can also point you to the entire OTU data set. so you could find high ID% matches with blastn.


## FILES INCLUDED: 
#alignment:3dom-org.18s28s.cmalign.f.fasta
#tree:3dom-org.18s28s.cmalign.f.fasta.contree
#tree building log with model parameters: 3dom-org.18s28s.cmalign.f.fasta.log
#taxonomic assigment table: tree604.notres.tassxstrings.archaea-bacteria-organelles-edit.tsv
