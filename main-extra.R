d.a66 <- 
combine_gene_sets_motif_targets(list(a66.set1, a66.set2, a66.set3, a66.set4, a66.set5, a66.set6), 20, "a66")
write.csv(d.a66[order(-score)], file = "../pip3-rna-seq-output/rds/wt-targets.csv")


d.a66.motif40 <- 	
combine_motif_targets(list(a66.set1, a66.set2, a66.set3, a66.set4, a66.set5, a66.set6), 40, "a66-no-facet")

plot_genes_by_motif(d, T, "test2")

plot_genes_by_motif(d.const.mot40[!target %in% d.a66.motif40[,target] & (set == 1 | set == 2)], T, "test5")


d.mut <- 
combine_gene_sets_motif_targets(list(set1, set2, set3, set4, set5, set6, set7, set8), 20, "const")
d.mut <- d.mut[!d.mut$target %in% d.a66$target]
write.csv(d.mut[order(-score)], file = "../pip3-rna-seq-output/rds/mut-targets.csv")


d.const.mot40 <- 
combine_motif_targets(list(set1, set2, set3, set4, set5, set6, set7, set8), 40, "const-no-facet")

d <- 
combine_motif_targets(list(egf.lrt), 40, "egf-no-facet")

