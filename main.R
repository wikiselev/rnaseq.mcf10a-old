# test of change remote

source("functions.R")

# PRELIMINARY CALCULATIONS process raw reads from htseq
process_raw_read_counts()
# perform differential expression analysis
differential_expression()

# ISMARA ANALYSIS process ISMARA report data
plot_ismara_averaged_activities()
plot_ismara_activities()
merge_ismara_target_files_new()
# Figure 4 in the paper
heatmap_ismara_activities()

# INITIALIZATION initialize gene sets of interest
initialize_sets()
# plot venn diagrams
create_venn_diagrams()

# ANALYSIS OF WT CELL LINES cluster mutational gene sets
clust_wt()
# get the clusters
clusts.wt <- get_clust_wt()
# GO analysis of mutational gene sets and their clusters
go_clust_wt(clusts.wt)
# plot GO plots for clusters
plot_go_clust_wt()
# TF binding motif differential activity analysis of WT gene sets
motif_diff_wt(clusts.wt)

# ANALYSIS OF MUTATIONAL CELL LINES cluster mutational gene sets
clust_mut()
# get the clusters
clusts.mut <- get_clust_mut()
# GO analysis of mutational gene sets and their clusters
go_clust_mut(clusts.mut)
# plot GO plots for clusters
plot_go_clust_mut()
# TF binding motif differential activity analysis of mutational gene sets
motif_diff_mut(clusts.mut)


# create data files for Vero
files <- list.files("../pip3-rna-seq-output/GO/")
for (f in files) {
    go_genes_for_vero(f)
}

# analysis of similarities of our results with the results of a new paper: Hart, J. R. et al. The
# butterfly effect in cancer: A single base mutation can remodel the cell.  Proc. Natl. Acad. Sci. U. S.
# A. 112, 1131–1136 (2015).
butterfly_paper_comparisons()

prepare_data_for_len_phil_effect()

# table 2 in the paper
d <- read.csv("../pip3-rna-seq-output/rds/prdm1-targets.csv")
plot_prdm1_genes(d$Motif.target.gene, T, "prdm1-targets")

export_raw_counts_for_submission() 
