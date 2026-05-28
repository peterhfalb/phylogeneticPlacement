## R script to prune 16S tree

library(ape)


tree <- read.tree("results/FAB2_2022_16S_phyloPlacement.newick") # replace with the name of your tree file

# Read query tip names directly from the filtered fasta
fasta_lines <- readLines("results/FAB2_2022_16S_phyloPlacement_aligned_masked_filtered.fasta") # replace with the name of your fasta file
query_tips <- sub("^>", "", fasta_lines[startsWith(fasta_lines, ">")])

# prune the tree
pruned_tree <- keep.tip(tree, query_tips)
write.tree(pruned_tree, "results/FAB2_2022_16S_bacteriaTree_final_unrooted.newick")

# root the tree and make it ultrametric
rooted_tree <- midpoint.root(pruned_tree) # root the tree at the midpoint first
pruned_tree_ultra <- chronos(rooted_tree)
write.tree(rooted_tree, "results/FAB2_2022_16S_bacteriaTree_final_rooted.newick")
write.tree(pruned_tree_ultra, "results/FAB2_2022_16S_bacteriaTree_final_ultrametric.newick")