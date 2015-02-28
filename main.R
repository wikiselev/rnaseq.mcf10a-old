source("functions.R")

##########################
# PRELIMINARY CALCULATIONS
#
# process raw reads from htseq
process_raw_read_counts()
# process mirna reads
mirna_read_process()
# import miRTarBase database with miRNA targets
miRTarBase()
# perform differential expression analysis
differential_expression()
##########################

################
# INITIALIZATION
#
# initialize gene sets of interest
initialize_sets()
# plot venn diagrams
create_venn_diagrams()
################

###################################
# ANALYSIS OF WT CELL LINES
#
# cluster mutational gene sets
clust_wt()
# get the clusters
clusts.wt <- get_clust_wt()
# GO analysis of mutational gene sets and their clusters
go_clust_wt(clusts.wt)
# plot GO plots for clusters
plot_go_clust_wt()
# plot GO trees for clusters
quickgo_clust_wt()
# TF binding motif differential activity analysis of WT gene sets
motif_diff_wt(clusts.wt)
###################################

###################################
# ANALYSIS OF MUTATIONAL CELL LINES
#
# make plots for general overview of gene expression at 0 min in mutated cell lines
genes <- c("FOXO1", "FOXO3", "FOXO4", "FOXO6", "JUN", "JUNB")
mut_overview_custom(genes)
mut_overview_sig()
# cluster mutational gene sets
clust_mut()
# get the clusters
clusts.mut <- get_clust_mut()
# GO analysis of mutational gene sets and their clusters
go_clust_mut(clusts.mut)
# plot GO plots for clusters
plot_go_clust_mut()
# plot GO trees for clusters
quickgo_clust_mut()
# TF binding motif differential activity analysis of mutational gene sets
motif_diff_mut(clusts.mut)
# plot prdm1 genes for publication
prdm1_genes_figure()
###################################

# respiration clusters
g.main83 <- get_genes_from_go("main-8-3", 10)
plot_genes(g.main83, T, "go-main-8-3")
g.main61 <- get_genes_from_go("main-6-1", 10)
plot_genes(g.main61, T, "go-main-6-1")
g.a6621 <- get_genes_from_go("a66-2-1", 10)
plot_genes(g.a6621, T, "go-a66-2-1")
venn(list("main-8-3" = g.main83, "main-6-1" = g.main61, "a66-2-1" = g.a6621), T, "respiration")

# translation cluster
g.a6623 <- get_genes_from_go("a66-2-3", 10)
plot_genes(g.a6623, T, "go-a66-2-3")

# splicing
g.a6622 <- get_genes_from_go("a66-2-2", 10)
plot_genes(g.a6622, T, "go-a66-2-2")

# ncRNA, rRNA, mRNA metabolism
g.a6611 <- get_genes_from_go("a66-1-1", 10)
plot_genes(g.a6611, T, "go-a66-1-2")
g.main71 <- get_genes_from_go("main-7-1", 7)
plot_genes(g.main71, T, "go-main-7-1")
g.main32 <- get_genes_from_go("main-3-2", 6)
plot_genes(g.main32, T, "go-main-3-2")
venn(list("main-7-1" = g.main71, "main-3-2" = g.main32, "a66-1-1" = g.a6611), T, "respiration")


plot_mirna(mirna.a66, "a66")
plot_mirna(mirna.egf, "egf")
plot_mirna(mirna.a66[mirna.a66 %in% mirna.egf], "a66-egf")
plot_mirna(mirna.a66[!mirna.a66 %in% mirna.egf], "a66-negf")

d <- get_miRTarBase()
mirna.targets <- unique(d[miRNA %in% mirna.a66, Target.Gene])

mirna.targets.ens <- hgnc_symbol_to_ensembl_id(mirna.targets)[,1]

venn(list("miRNA targets" = mirna.targets.ens, "set1" = set1, "set2" = set2), T, "mirna-set1-set2")
venn(list("miRNA targets" = mirna.targets.ens, "a66-1" = a66.set1, "a66-2" = a66.set2), T, "mirna-a66-1-a66-2")
venn(list("miRNA targets" = mirna.targets.ens, "EGF responding" = egf.lrt), T, "mirna-egf")

plot_genes(mirna.targets.ens[mirna.targets.ens %in% set1], F, "mirna-tar-set1")
plot_genes(mirna.targets.ens[mirna.targets.ens %in% set2], F, "mirna-tar-set2")

plot_genes(mirna.targets.ens[mirna.targets.ens %in% a66.set1], F, "mirna-tar-a66-1")
plot_genes(mirna.targets.ens[mirna.targets.ens %in% a66.set2], F, "mirna-tar-a66-2")

# d <- motif_diff_activity(mirna.targets.ens[mirna.targets.ens %in% a66.set1],
# 	mirna.targets.ens[mirna.targets.ens %in% a66.set2], "mirna-tar-a66-1", "mirna-tar-a66-2", "mirna-a66")

GO(mirna.targets.ens[mirna.targets.ens %in% set1], all.genes, "mirna-tar-set1", 0.05)
GO(mirna.targets.ens[mirna.targets.ens %in% set2], all.genes, "mirna-tar-set2", 0.05)
GO(mirna.targets.ens[mirna.targets.ens %in% a66.set1], all.genes, "mirna-tar-a66-1", 0.05)
GO(mirna.targets.ens[mirna.targets.ens %in% a66.set2], all.genes, "mirna-tar-a66-2", 0.05)
GO(mirna.targets.ens, all.genes, "mirna-tar", 0.05)

mirna.set1 <-unique(d[Target.Gene %in% ensembl_id_to_hgnc_symbol(mirna.targets.ens[mirna.targets.ens %in% set1])[,2], miRNA])
mirna.set1[mirna.set1 %in% mirna.a66]
"hsa-miR-7-5p"   "hsa-miR-197-3p" "hsa-miR-27b-5p"

# with clust_boot_a66 I checked if clustering would be different if I only 
# consider a66 and wt time courses excluding pten and ki
# It was more or less the same - so I am not doing it in the final version of
# the manuscript

# clust_boot_a66(a66_1[a66_1 %in% egf.lrt], 2, 6, "a66_only", "1")
# clust_boot_a66(a66_1[!a66_1 %in% egf.lrt], 2, 6, "a66_only", "2")
# clust_boot_a66(a66_2[a66_2 %in% egf.lrt], 2, 6, "a66_only", "3")
# clust_boot_a66(a66_2[!a66_2 %in% egf.lrt], 2, 6, "a66_only", "4")
# clust_boot_a66(a66_3[a66_3 %in% egf.lrt], 2, 6, "a66_only", "5")
# clust_boot_a66(a66_3[!a66_3 %in% egf.lrt], 2, 6, "a66_only", "6")

# plot_bootstrap_data("a66_only")

# clust1 <- get_clust_genes("a66_only", "1-2")
# clust2 <- get_clust_genes("a66_only", "2-2")
# clust3 <- get_clust_genes("a66_only", "3-2")
# clust4 <- get_clust_genes("a66_only", "4-3")
# clust5 <- get_clust_genes("a66_only", "5-2")
# clust6 <- get_clust_genes("a66_only", "6-3")


# clusts <- list(clust1$partition, clust2$partition, clust3$partition,
# 	clust4$partition, clust5$partition, clust6$partition)
# inds <- c(1, 2, 3, 4, 5, 6)

# plot_all_clusts_a66(clusts, "a66_only", inds)

# ind <- 1
# for(j in clusts) {
# 	for (i in 1:length(unique(j))) {
# 		GO(names(j[j == i]), all.genes, paste0("a66-", inds[ind], "-", i), 0.01)
# 	}
# 	# print(inds[ind])
# 	ind <- ind + 1
# }







# search for semantic similarities
# go11 <- c("GO:0042221", "GO:0006796", "GO:0006468", "GO:0006629", "GO:0042311", "GO:0008015", "GO:0050896", "GO:0016310", "GO:0051923", "GO:0035284")
# go12 <- c("GO:0006928", "GO:0042221", "GO:0072358", "GO:0072359", "GO:0040011", "GO:0009653", "GO:0001944", "GO:0045597", "GO:0060396", "GO:0001568")
# go13 <- c("GO:0042060", "GO:0010631", "GO:0007596", "GO:0043542", "GO:0007599", "GO:0009611", "GO:0050974", "GO:0010757", "GO:0045785", "GO:0030168")

# mgoSim(go11, go12, ont = "BP", measure = "Wang", combine = "BMA")

# plot all interesting genes with constitutive mutations effect sorted by the
# effect significance -- this had been done to find out 10-20 the most interesting
# genes and do a qPCR validation on them
filtering_of_constitutive_effects(0.05)

# genes selected by me, Vero and Nicolas based on filtering of const effects
int.genes <- c("ACPL2", "ATP9A", "AXL", "BCL7A", "CLASP1", "GLG1", "GPC1", "KRT10", "LPAR3", "NTN4", "RGCC", "SDF2L1", "THBD", "CCAT1", "DAW1", "PI3", "ABCC2", "TACC2", "TRIML2", "NRG4", "TNNI2", "ZNF584", "IL27RA", "EPGN", "STK40", "PIK3R1", "PROM2", "BAK1", "VAV2", "CAPN6", "VAV1", "KRT13", "MYEOV", "PLA2G4F", "CES1", "AKR1B15", "PTK6", "S100A8", "PREX1", "H19", "S100A9", "SYT8", "NQO1", "HMGN5", "MED24", "ATP5G3", "NGFRAP1", "DTX3", "LBH", "ATP8B2", "ID4", "COL4A3", "CD27-AS1", "CALHM2", "MYL9", "ASPA", "VWA5A", "PXDN", "TMEM139", "STON1", "FRMD3", "LMCD1", "CLDN1", "COL4A4", "LPXN", "AKAP12", "F2R", "PIK3IP1", "ARHGAP29", "CAMK2N1", "INPP4B", "LIF", "SLC46A3", "KIF3C", "DENND2D", "NNMT", "ARHGEF28", "REPS2", "HPSE", "SLFN5", "HEG1", "CREB3L4", "C16orf62", "TCP11L1", "PTK2", "RASSF5", "RNF152", "TPM1", "HIP1", "CALCOCO1", "SWAP70", "CAMLG", "SPTBN1", "PDLIM1", "HERPUD2", "RGS11", "SPARC", "CDH13", "SRPX", "C8orf47", "QPCT", "NME5", "MAGEH1", "SLIT3", "ZNF37", "BDH1", "THNSL2", "MCAM", "PLA2R1", "ZNF585A", "CFI", "ZNF358", "FN1", "ITGB2", "NMNAT2", "NMNAT3", "DOCK10", "DAB2", "ATL1", "TUBA1A", "LOX", "ARG2", "CARD16", "SCEL", "FAP", "FSTL1", "LTBP1", "PKN1", "DIXDC1", "CRIM1", "AKT3", "APLP2")
# genes proposed for further pcr experiment
pcr.genes <- c("DENND2D", "ID4", "INPP4B", "LPXN", "PIK3IP1", "PTK2", "TMEM139", "VWA5A", "CCAT1", "PIK3R1", "PREX1", "PROM2", "VAV2", "AXL", "GPC1", "LPAR3", "AKT3", "FN1", "NMNAT2", "NMNAT3")

plot_genes(pcr.genes, F, "pcr")
plot_genes(pcr.genes, T, "pcr")

# genes selected by Len based on filtering of const effects and proposed for
# further pcr experiment
# all of these genes are actually in int.genes!!!
len.grp1 <- c("ARHGEF28", "RGCC", "PTK2", "CALCOCO1", "LBH", "PIK3IP1", "HERPUD2")
len.grp2 <- c("STON1", "HIP1")
len.grp3 <- c("STK40", "EPGN", "DAW1", "BAK1")
len.grp4 <- c("INPP4B", "CLDN1")

plot_genes(len.grp1, F, "len-grp1")
plot_genes(len.grp1, T, "len-grp1")
plot_genes(len.grp2, F, "len-grp2")
plot_genes(len.grp2, T, "len-grp2")
plot_genes(len.grp3, F, "len-grp3")
plot_genes(len.grp3, T, "len-grp3")
plot_genes(len.grp4, F, "len-grp4")
plot_genes(len.grp4, T, "len-grp4")

# venn diagram comparing genes selected by us and by Len
venn(list("VVN" = pcr.genes, "LP" = c(len.grp1, len.grp2, len.grp3, len.grp4)), T, "pcr-genes")
pcr.genes[pcr.genes %in% c(len.grp1, len.grp2, len.grp3, len.grp4)]
# "INPP4B"  "PIK3IP1" "PTK2"

# ISMARA ANALYSIS
#
# process ISMARA report data
plot_ismara_averaged_activities()
plot_ismara_activities()
merge_ismara_target_files()

# analysis of motif differential activity using Kolmogorov-Smirnov and
# Wilcox tests
d <- motif_diff_activity(set1, set2, "dynamic", "static", "main")
plot_motifs_tfs("SRF.p3")
plot_motifs_tfs("PRDM1.p3")

# svd analysis of motif activities obtained from ismara report.
# need to add reordering to activity matrix!!! otherwise get wrong results!
# reordering should be done by arrange_count_matrix_columns function, but
# somehow it does not work on ISMARA activity table - wt columns are messed up
d <- s_v_d_ismara()

plot_motifs_tfs("PRDM1.p3")
t <- get_ismara_targets("PRDM1.p3")
g <- unique(t[,target_ens])
g <- g[g %in% ki & g %in% pten]
plot_genes(g, F, "PRDM1_targets")

GO(g, "ismara_PRDM1.p3_targets")

int.motifs <- c("PRDM1.p3", "HBP1_HMGB_SSRP1_UBTF.p2", "bHLH_family.p2",
	"FOX.D1.D2..p2", "ZNF384.p2", "AHR_ARNT_ARNT2.p2", "UUUUUGC")
plot_motifs_tfs("PRDM1.p3")
plot_motifs_tfs("HBP1_HMGB_SSRP1_UBTF.p2")
plot_motifs_tfs("bHLH_family.p2")
plot_motifs_tfs("FOX.D1.D2..p2")
plot_motifs_tfs("ZNF384.p2")
plot_motifs_tfs("AHR_ARNT_ARNT2.p2")
plot_motifs_tfs("UUUUUGC")

t <- get_ismara_targets(int.motifs)
g <- unique(t[,target_ens])

####
# time course analysis
order_time_courses_by_wt_peaks()
plot_peaks_old()
####

plot_motif_diff_targets(d, "PRDM1.p3")

d <- plot_motif_target_hist(ki[!(ki %in% pten)], 5)
d$p
d <- plot_motif_target_hist(pten[!(pten %in% ki)], 5)
d$p

split_targets_by_motifs(a66, 1, 5)


# alternative splicing
# this should be done on a cluster because it takes several hours to finish
# please use script e.g. "Rscript splicing_on_cluster.R pten" to compare wt and pten
alt_splicing("wt", "ki", 0, 0)
alt_splicing("wt", "pten", 0, 0)

# find out which of the genes are lincRNAs
linc_rna("mut")
linc_rna("a66")
# plot genes selected by Vero which might be regulated by lincRNAs (affected by mutations)
# need to rerun it on new lists - in the old ones some gene names corresponded to
# hundreds of ENSEMBL IDs!!!
plot_linc_rna()

# create data files for Vero
files <- list.files("../pip3-rna-seq-output/GO/")
for(f in files) {
	go_genes_for_vero(f)
}

# analysis of similarities of our results with the results of a new
# paper: Hart, J. R. et al. The butterfly effect in cancer:
# A single base mutation can remodel the cell.
# Proc. Natl. Acad. Sci. U. S. A. 112, 1131–1136 (2015).
butterfly_paper_comparisons()
heatmap_ismara_activities()
prepare_data_for_len_phil_effect()
