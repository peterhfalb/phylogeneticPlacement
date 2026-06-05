## R script to prune 18S eukaryote placement tree

library(ape)


tree <- read.tree("results/nameOfOutputFiles.newick") # replace with the name of your tree file

# Read query tip names directly from the filtered fasta
fasta_lines <- readLines("results/nameOfOutputFiles_aligned_padded_filtered.fasta") # replace with the name of your fasta file
query_tips <- sub("^>", "", fasta_lines[startsWith(fasta_lines, ">")])

# prune the tree
pruned_tree <- keep.tip(tree, query_tips)
write.tree(pruned_tree, "results/nameOfOutputFiles_eukaryoteTree_final_unrooted.newick")

# root the tree and make it ultrametric
rooted_tree <- midpoint.root(pruned_tree) # root the tree at the midpoint first
pruned_tree_ultra <- chronos(rooted_tree)
write.tree(rooted_tree, "results/nameOfOutputFiles_eukaryoteTree_final_rooted.newick")
write.tree(pruned_tree_ultra, "results/nameOfOutputFiles_eukaryoteTree_final_ultrametric.newick")
