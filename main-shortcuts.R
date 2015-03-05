process_raw_read_counts <- function() {
  # create a count matrix from raw counts produced by htseq_count
  files <- list.files("../pip3-rna-seq-input/htseq-read-count-raw", full.names = FALSE)
  count.matrix <- count_matrix_from_files(files, F)
  count.matrix <- arrange_count_matrix_columns(count.matrix)
  # plot correlation heatmap of count.matrix
  cor_heatmap(count.matrix, "cor-mat-raw-with-bad-rep.pdf")

  # remove bad replicates: wt_3_40 and wt_4_40
  files <- files[!grepl("wt_3_40", files)]
  files <- files[!grepl("wt_4_40", files)]
  files <- files[!grepl("wt_3_0", files)]

  count.matrix <- count_matrix_from_files(files, F)
  count.matrix <- arrange_count_matrix_columns(count.matrix)
  # plot correlation heatmap of count.matrix
  cor_heatmap(count.matrix, "cor-mat-raw.pdf")
  # now correlation heatmap looks much better!

  # change colnames for b replicates from "wtb" to "wt"
  cn <- colnames(count.matrix)
  cn[grepl("wtb", cn)] <- 
    c("wt_3b_0", "wt_3b_40", "wt_4b_40")
  colnames(count.matrix) <- cn

  # I talked to Simon Andrews - he said that removing lowly-expressed
  # genes by rowsums or any other method can be risky and DESeq2 will anyway
  # remove all those genes using negative binomial distribution.
  # So, finally I only remove genes which are not expressed in any condition
  # (count.matrix[rowSums(count.matrix) > 0, ]) and this will be the data set
  # which I use in the further analysis.

  # remove all 0-expressed genes
  count.matrix <- count.matrix[rowSums(count.matrix) > 0, ]

  # convert count.matrix into data.table
  dat <- count_matrix_2_data_table(count.matrix)
  # check read density distribution and cumulative distribution
  plot_read_density(dat, "raw")
  plot_read_ecdf(dat, "raw")
  # read distributions of different samples are different
  # let's now use a count.matrix normalized by the library size

  # normalize read counts by library size using DESeq2 estimateSizeFactors function
  count.matrix <- count_matrix_from_files(files, T)
  count.matrix <- arrange_count_matrix_columns(count.matrix)
  # plot correlation heatmap of count.matrix
  cor_heatmap(count.matrix, "cor-mat-raw-norm.pdf")
  # normalization does not change correlation -- this is a sanity check!

  # change colnames for b replicates from "wtb" to "wt"
  cn <- colnames(count.matrix)
  cn[grepl("wtb", cn)] <- 
    c("wt_3b_0", "wt_3b_40", "wt_4b_40")
  colnames(count.matrix) <- cn
  # remove all 0-expressed genes
  count.matrix <- count.matrix[rowSums(count.matrix) > 0, ]
  # convert count.matrix into data.table
  dat <- count_matrix_2_data_table(count.matrix)
  # check read density distribution and cumulative distribution
  plot_read_density(dat, "norm")
  plot_read_ecdf(dat, "norm")
  # now distributions are much better! Normalized count.matrix can be saved:
  saveRDS(count.matrix, "../pip3-rna-seq-output/rds/count-matrix.rds")
  # scale and center count matrix -- needed for svd analysis and for clustering
  count.matrix.scaled <- t(scale(t(count.matrix), center = T, scale = T))
  saveRDS(count.matrix.scaled, "../pip3-rna-seq-output/rds/count-matrix-scaled.rds")

  # principal component analysis
  pca()

  # annotate genes in the count.matrix with associated gene names
  # file gene_names_GRCh37.p13.txt was downloaded from Ensembl Biomart on 06/07/14
  gene.ann <- read.csv("../pip3-rna-seq-input/annotations/gene-names-GRCh37.p13.txt")
  colnames(gene.ann) <- c("ensembl_gene_id", "hgnc_symbol")
  res <- data.frame(ensembl_gene_id = rownames(count.matrix))
  res <- merge(res, gene.ann)
  saveRDS(res, "../pip3-rna-seq-output/rds/count-matrix-ann.rds")
  
  # I also average all replicates at each condition and each time point
  ave_repl(count.matrix, "")
  ave_repl(count.matrix.scaled, "-scaled")

  # normalize genes by their lengths
  # file gene_lengths_GRCh37.p13.txt was downloaded from Ensembl Biomart on 06/07/14
  gene.len <- read.csv("../pip3-rna-seq-input/annotations/gene-lengths-GRCh37.p13.txt")
  colnames(gene.len) <- c("ensembl_gene_id", "gene_start", "gene_end")
  gene.len$gene_length <- abs(gene.len$gene_end - gene.len$gene_start)
  res <- data.frame(ensembl_gene_id = rownames(count.matrix))
  res <- merge(res, gene.len[ , c(1, 4)])

  # I multiply everything by 100 to shift distribution of values to 0:
  # hist(log10( count.matrix / res$gene_length * 100))
  count.matrix.by.length <- count.matrix / res$gene_length * 100
  res <- average_count_matrix(count.matrix.by.length, conditions)
  plot.data.all <- create_plot_table(res[[1]], res[[2]])
  saveRDS(plot.data.all, "../pip3-rna-seq-output/rds/plot-time-courses-all-genes-norm-by-length.rds")
}

initialize_sets <- function() {
	# define initial gene sets
	egf.lrt <<- rownames(get_diff_expr("wt", 0.01))

	egf.wt15 <<- rownames(get_diff_expr("wt-0-wt-15", 0.05))
	egf.wt40 <<- rownames(get_diff_expr("wt-0-wt-40", 0.05))
	egf.wt90 <<- rownames(get_diff_expr("wt-0-wt-90", 0.05))
	egf.wt180 <<- rownames(get_diff_expr("wt-0-wt-180", 0.05))
	egf.wt300 <<- rownames(get_diff_expr("wt-0-wt-300", 0.05))

	ki <<- rownames(get_diff_expr("wt-0-ki-0", 0.05))
	pten <<- rownames(get_diff_expr("wt-0-pten-0", 0.05))
	a66 <<- rownames(get_diff_expr("wt-0-ko-0", 0.05))
	pi3k.pten <<- rownames(get_diff_expr("pten-0-ki-0", 0.05))

	egf.pten <<- rownames(get_diff_expr("pten", 0.01))
	egf.ki <<- rownames(get_diff_expr("ki", 0.01))
	egf.a66 <<- rownames(get_diff_expr("ko", 0.01))

	pten.tc <<- rownames(get_diff_expr("wt-pten", 0.01))
	ki.tc <<- rownames(get_diff_expr("wt-ki", 0.01))
	a66.tc <<- rownames(get_diff_expr("wt-ko", 0.01))

	pten.tc.cond <<- rownames(get_diff_expr("wt-pten-cond", 0.01))
	ki.tc.cond <<- rownames(get_diff_expr("wt-ki-cond", 0.01))
	a66.tc.cond <<- rownames(get_diff_expr("wt-ko-cond", 0.01))

	pi3k.related <<- unique(c(pten.tc.cond, ki.tc.cond))

	a66_nost_eff <<- rownames(get_diff_expr("ko-0-konost-300", 0.05))

	mirna.egf <<- rownames(get_diff_expr("mirna-wt", 0.01))
	mirna.a66 <<- rownames(get_diff_expr("mirna-ko-0-konost-300", 0.05))
	# length(a66_nost_eff)
	# 4168 - a lot of genes were affected by A66 treatment in the absence of EGF

	# new gene universe!!!
	all.genes <<- unique(c(egf.lrt, ki, pten, pi3k.related, a66.tc.cond, a66_nost_eff))

	# constitutive mutation gene sets

	# 711 genes
	set1 <<- ki[ki %in% pten & ki %in% egf.lrt & ki %in% pi3k.related]
	# 407 genes
	set2 <<- ki[ki %in% pten & !ki %in% egf.lrt & ki %in% pi3k.related]
	# 802 genes
	set3 <<- ki[!ki %in% pten & ki %in% egf.lrt & ki %in% pi3k.related]
	# 416 genes
	set4 <<- ki[!ki %in% pten & !ki %in% egf.lrt & ki %in% pi3k.related]
	# 1384 genes
	set5 <<- pten[!pten %in% ki & pten %in% egf.lrt & pten %in% pi3k.related]
	# 1095 genes
	set6 <<- pten[!pten %in% ki & !pten %in% egf.lrt & pten %in% pi3k.related]
	# 3183 genes
	set7 <<- egf.lrt[!egf.lrt %in% pten & !egf.lrt %in% ki & egf.lrt %in% pi3k.related]
	# 3002 genes
	set8 <<- pi3k.related[!pi3k.related %in% pten & !pi3k.related %in% ki & !pi3k.related %in% egf.lrt]
	# 2493 genes
	set9 <<- egf.lrt[!egf.lrt %in% pten & !egf.lrt %in% ki & !egf.lrt %in% pi3k.related]

	# a66 gene sets

	a66.set1 <<- a66_nost_eff[!a66_nost_eff %in% a66.tc.cond & a66_nost_eff %in% egf.lrt]
	a66.set2 <<- a66_nost_eff[!a66_nost_eff %in% a66.tc.cond & !a66_nost_eff %in% egf.lrt]
	a66.set3 <<- a66.tc.cond[!a66.tc.cond %in% a66_nost_eff & a66.tc.cond %in% egf.lrt]
	a66.set4 <<- a66.tc.cond[!a66.tc.cond %in% a66_nost_eff & !a66.tc.cond %in% egf.lrt]
	a66.set5 <<- a66_nost_eff[a66_nost_eff %in% a66.tc.cond & a66_nost_eff %in% egf.lrt]
	a66.set6 <<- a66_nost_eff[a66_nost_eff %in% a66.tc.cond & !a66_nost_eff %in% egf.lrt]
}

differential_expression <- function() {
	# find EGF responding genes in all conditions
	diff_expr_time_course("wt")
	diff_expr_time_course("pten")
	diff_expr_time_course("ki")
	diff_expr_time_course("ko")
	# calculate differentially expressed genes between all conditions and WT at
	# 0 time point
	diff_expr_pairwise("wt", "ko", 0, 0)
	diff_expr_pairwise("wt", "ki", 0, 0)
	diff_expr_pairwise("wt", "pten", 0, 0)
	diff_expr_pairwise("pten", "ki", 0, 0)
	# calculate differentially expressed genes using A66 control 
	# (no EGF stimulation) at 300 min
	diff_expr_pairwise("ko", "konost", 0, 300)

	# check interactions between all conditions and wt
	diff_expr_two_time_courses("wt", "pten")
	diff_expr_two_time_courses("wt", "ki")
	diff_expr_two_time_courses("wt", "ko")

	# check condition+interaction between all conditions and wt
	diff_expr_two_time_courses_cond("wt", "pten")
	diff_expr_two_time_courses_cond("wt", "ki")
	diff_expr_two_time_courses_cond("wt", "ko")

	# mirna differential expression
	diff_expr_time_course_mirna("wt")
	diff_expr_pairwise_mirna("ko", "konost", 0, 300)
}

create_venn_diagrams <- function() {
	venn(list(EGF = egf.lrt, "PI3K KI" = ki, "PTEN KO" = pten), T,
		"const-mut-vs-egf-lrt")
	venn(list(WT = egf.lrt, A66 = egf.a66, "PTEN KO" = egf.pten, "PI3K KI" = egf.ki),
		T, "all-egf")
	venn(list(WT = egf.lrt, A66 = egf.a66, "PTEN KO" = egf.pten, "PI3K KI" = egf.ki),
		F, "all-egf1")
	venn_ellipses(list("EGF responding in WT" = egf.lrt,
			  "PTEN KO + PI3K KI, PI3K related" = pi3k.related,
			  "PI3K KI, constitutive" = ki,
			  "PTEN KO, constitutive" = pten), F, "const-mut-tc-vs-egf-lrt")
	venn(list("A66 effect (no EGF)" = a66_nost_eff,
			  "PTEN KO + PI3K KI, PI3K related" = pi3k.related,
			  "PI3K KI, constitutive" = ki,
			  "PTEN KO, constitutive" = pten), F, "const-mut-tc-vs-a66-no-egf")

	venn(list("A66 effect (no EGF)" = a66_nost_eff,
			  "A66 effect" = a66.tc.cond, "EGF responding in WT" = egf.lrt),
			  T, "a66-egf-tc-vs-a66-no-egf")
}

mut_overview_custom <- function(genes){
	# import averaged (by replicates) and normalized count matrix
	cm <- readRDS("../pip3-rna-seq-output/rds/count-matrix-av.rds")
	# select WT and both mutations at 0 time point
	cm <- cm[,c(1,2,14,20)]
	cm$gname <- ensembl_id_to_hgnc_symbol(cm$id)$hgnc_symbol
	cm <- as.data.frame(cm)

	cm$wt_0 <- log2(cm$wt_0 + 1)
	cm$pten_0 <- log2(cm$pten_0 + 1)
	cm$ki_0 <- log2(cm$ki_0 + 1)

	p <- ggplot(cm, aes(wt_0, pten_0)) + 
	  geom_point(colour="lightblue") +
	  geom_point(data = cm[cm$gname %in% genes, ], aes(x=wt_0, y=pten_0), colour="red", size=4) +
	  geom_text(data = cm[cm$gname %in% genes, ], aes(label = gname), hjust=-0.2, vjust=-0.7, size = 4) +
	  geom_abline(intercept = 0, slope = 1, colour = "black", size = 1) +
	  theme_bw()
	pdf(file = "../pip3-rna-seq-output/figures/mut-0-overview-pten.pdf", w = 8, h = 8)
	print(p)
	dev.off()

	p <- ggplot(cm, aes(wt_0, ki_0)) + 
	  geom_point(colour="lightblue") +
	  geom_point(data = cm[cm$gname %in% genes, ], aes(x=wt_0, y=ki_0), colour="red", size=4) +
	  geom_text(data = cm[cm$gname %in% genes, ], aes(label = gname), hjust=-0.2, vjust=-0.7, size = 4) +
	  geom_abline(intercept = 0, slope = 1, colour = "black", size = 1) +
	  theme_bw()
	pdf(file = "../pip3-rna-seq-output/figures/mut-0-overview-ki.pdf", w = 8, h = 8)
	print(p)
	dev.off()

	p <- ggplot(cm, aes(pten_0, ki_0)) + 
	  geom_point(colour="lightblue") +
	  geom_point(data = cm[cm$gname %in% genes, ], aes(x=pten_0, y=ki_0), colour="red", size=4) +
	  geom_text(data = cm[cm$gname %in% genes, ], aes(label = gname), hjust=-0.2, vjust=-0.7, size = 4) +
	  geom_abline(intercept = 0, slope = 1, colour = "black", size = 1) +
	  theme_bw()
	pdf(file = "../pip3-rna-seq-output/figures/mut-0-overview-pten-ki.pdf", w = 8, h = 8)
	print(p)
	dev.off()
}

mut_overview_sig <- function(){
	# import averaged (by replicates) and normalized count matrix
	cm <- readRDS("../pip3-rna-seq-output/rds/count-matrix-av.rds")
	# select WT and both mutations at 0 time point
	cm <- cm[,c(1,2,14,20)]
	# cm$gname <- ensembl_id_to_hgnc_symbol(cm$id)$hgnc_symbol
	cm <- as.data.frame(cm)

	cm$wt_0 <- log2(cm$wt_0 + 1)
	cm$pten_0 <- log2(cm$pten_0 + 1)
	cm$ki_0 <- log2(cm$ki_0 + 1)

	p <- ggplot(cm, aes(wt_0, pten_0)) + 
	  geom_point(colour="gray") +
	  geom_point(data = cm[cm$id %in% pten, ], aes(x=wt_0, y=pten_0), colour="blue") +
	  geom_abline(intercept = 0, slope = 1, colour = "black", size = 1) +
	  theme_bw()
	ggsave("../pip3-rna-seq-output/figures/mut-0-overview-pten.png", p, width=6, height=6)

	p <- ggplot(cm, aes(wt_0, ki_0)) + 
	  geom_point(colour="gray") +
	  geom_point(data = cm[cm$id %in% ki, ], aes(x=wt_0, y=ki_0), colour="green") +
	  geom_abline(intercept = 0, slope = 1, colour = "black", size = 1) +
	  theme_bw()
	ggsave("../pip3-rna-seq-output/figures/mut-0-overview-ki.png", p, width=6, height=6)

	p <- ggplot(cm, aes(ki_0, pten_0)) + 
	  geom_point(colour="gray") +
	  geom_point(data = cm[cm$id %in% pi3k.pten, ], aes(x=ki_0, y=pten_0), colour="black") +
	  # geom_point(data = cm[cm$id %in% ki, ], aes(x=ki_0, y=pten_0), colour="green", size = 1.5) +
	  geom_abline(intercept = 0, slope = 1, colour = "black", size = 1) +
	  theme_bw()
	ggsave("../pip3-rna-seq-output/figures/mut-0-overview-pten-ki.png", p, width=6, height=6)
}

clust_mut <- function() {
	clust_boot(set1, 2, 8, "main", 1)
	clust_boot(set3, 2, 8, "main", 3)
	clust_boot(set5, 2, 8, "main", 5)
	clust_boot(set7, 2, 8, "main", 7)

	clust_boot(set2, 2, 6, "main", 2)
	clust_boot(set4, 2, 6, "main", 4)
	clust_boot(set6, 2, 6, "main", 6)
	clust_boot(set8, 2, 6, "main", 8)

	plot_bootstrap_data("main")
}

clust_wt <- function() {
	clust_boot(a66.set1, 2, 6, "a66", 1)
	clust_boot(a66.set2, 2, 6, "a66", 2)
	clust_boot(a66.set3, 2, 6, "a66", 3)
	clust_boot(a66.set4, 2, 6, "a66", 4)
	clust_boot(a66.set5, 2, 6, "a66", 5)
	clust_boot(a66.set6, 2, 6, "a66", 6)

	plot_bootstrap_data("a66")
}

get_clust_mut <- function() {
	clust1 <- get_clust_genes("main", "1-5")
	clust3 <- get_clust_genes("main", "3-3")
	clust5 <- get_clust_genes("main", "5-5")
	clust7 <- get_clust_genes("main", "7-4")

	clust2 <- get_clust_genes("main", "2-4")
	clust4 <- get_clust_genes("main", "4-2")
	clust6 <- get_clust_genes("main", "6-2")
	clust8 <- get_clust_genes("main", "8-3")

	clusts <- list(clust1$partition, clust2$partition, clust3$partition,
		clust4$partition, clust5$partition, clust6$partition, clust7$partition,
		clust8$partition)

	plot_all_clusts(clusts, "main", c(1, 2, 3, 4, 5, 6, 7, 8))
	return(clusts)
}

get_clust_wt <- function() {
	clust1 <- get_clust_genes("a66", "1-2")
	clust2 <- get_clust_genes("a66", "2-4")
	clust3 <- get_clust_genes("a66", "3-2")
	clust4 <- get_clust_genes("a66", "4-2")
	clust5 <- get_clust_genes("a66", "5-2")
	clust6 <- get_clust_genes("a66", "6-2")


	clusts <- list(clust1$partition, clust2$partition, clust3$partition,
		clust4$partition, clust5$partition, clust6$partition)

	plot_all_clusts(clusts, "a66", c(1, 2, 3, 4, 5, 6))
	return(clusts)
}

go_clust_mut <- function(clusts) {
	GO(set1, all.genes, "main-1", 0.05)
	GO(set2, all.genes, "main-2", 0.05)
	GO(set3, all.genes, "main-3", 0.05)
	GO(set4, all.genes, "main-4", 0.05)
	GO(set5, all.genes, "main-5", 0.05)
	GO(set6, all.genes, "main-6", 0.05)
	GO(set7, all.genes, "main-7", 0.05)
	GO(set8, all.genes, "main-8", 0.05)

	ind <- 1
	for(j in clusts) {
		for (i in 1:length(unique(j))) {
			GO(names(j[j == i]), all.genes, paste0("main-", c(1, 2, 3, 4, 5, 6, 7, 8)[ind], "-", i), 0.05)
		}
		# print(inds[ind])
		ind <- ind + 1
	}
}

go_clust_wt <- function(clusts) {
	GO(a66.set1, all.genes, "a66-1", 0.05)
	GO(a66.set2, all.genes, "a66-2", 0.05)
	GO(a66.set3, all.genes, "a66-3", 0.05)
	GO(a66.set4, all.genes, "a66-4", 0.05)
	GO(a66.set5, all.genes, "a66-5", 0.05)
	GO(a66.set6, all.genes, "a66-6", 0.05)

	ind <- 1
	for(j in clusts) {
		for (i in 1:length(unique(j))) {
			GO(names(j[j == i]), all.genes, paste0("a66-", c(1, 2, 3, 4, 5, 6)[ind], "-", i), 0.05)
		}
		# print(inds[ind])
		ind <- ind + 1
	}
}

motif_diff_mut <- function(clusts) {
	d.const12 <<- motif_diff_activity_new(set1, set2, "1", "2", "const-1-vs-2-new")
	d.const34 <<- motif_diff_activity_new(set3, set4, "3", "4", "const-3-vs-4-new")
	d.const56 <<- motif_diff_activity_new(set5, set6, "5", "6", "const-5-vs-6-new")
	d.const78 <<- motif_diff_activity_new(set7, set8, "7", "8", "const-7-vs-8-new")
	d.const13 <<- motif_diff_activity_new(set1, set3, "1", "3", "const-1-vs-3-new")
	d.const15 <<- motif_diff_activity_new(set1, set5, "1", "5", "const-1-vs-5-new")

	d.const12 <<- motif_diff_activity_new(set1, set2, "1", "2", "const-1-vs-2")
	d.const34 <<- motif_diff_activity_new(set3, set4, "3", "4", "const-3-vs-4")
	d.const56 <<- motif_diff_activity_new(set5, set6, "5", "6", "const-5-vs-6")
	d.const78 <<- motif_diff_activity_new(set7, set8, "7", "8", "const-7-vs-8")
	d.const13 <<- motif_diff_activity_new(set1, set3, "1", "3", "const-1-vs-3")
	d.const15 <<- motif_diff_activity_new(set1, set5, "1", "5", "const-1-vs-5")

	d <<- motif_diff_activity_new(names(clusts[[8]][clusts[[8]] == 1]),
		names(clusts[[8]][clusts[[8]] == 3]), "8-1", "8-3", "8-1-8-3")

	d <<- motif_diff_activity_new(names(clusts[[1]][clusts[[1]] == 2]),
		names(clusts[[2]][clusts[[2]] == 2]), "1-2", "2-2", "1-2-2-2")
}

motif_diff_wt <- function(clusts) {
	d.a6612 <<- motif_diff_activity_new(a66.set1, a66.set2, "1", "2", "a66-1-2")
	d.a6613 <<- motif_diff_activity_new(a66.set1, a66.set3, "1", "3", "a66-1-3")
	d.a6623 <<- motif_diff_activity_new(a66.set2, a66.set3, "2", "3", "a66-2-3")

	d.a661112 <<- motif_diff_activity_new(names(clusts[[1]][clusts[[1]] == 1]),
		names(clusts[[1]][clusts[[1]] == 2]), "1-1", "1-2", "a66-1-1-1-2")

	d <<- motif_diff_activity_new(names(clusts[[2]][clusts[[2]] == 1]),
		names(clusts[[2]][clusts[[2]] == 2]), "2-1", "2-2", "a66-2-1-2-2")

	d.a663132 <<- motif_diff_activity_new(names(clusts[[3]][clusts[[3]] == 1]),
		names(clusts[[3]][clusts[[3]] == 2]), "3-1", "3-2", "a66-3-1-3-2")

	d.a665152 <<- motif_diff_activity_new(names(clusts[[5]][clusts[[5]] == 1]),
		names(clusts[[5]][clusts[[5]] == 2]), "5-1", "5-2", "a66-5-1-5-2")

	d <<- motif_diff_activity_new(names(clusts[[1]][clusts[[1]] == 1]),
		names(clusts[[3]][clusts[[3]] == 2]), "1-1", "3-2", "a66-1-1-3-2")

	d <<- motif_diff_activity_new(names(clusts[[1]][clusts[[1]] == 2]),
		names(clusts[[3]][clusts[[3]] == 1]), "1-2", "3-1", "a66-1-2-3-1")
}

prdm1_genes_figure <- function() {
	d <- d.const12[motif == "PRDM1.p3" & diff==1, list(target, score)]
	d <- d[order(-score)]
	plot_prmd1_genes(d[,target], T, "prdm1-genes-figure")
}

plot_go_clust_mut <- function() {
	main_1 <<- annotate_go_terms("main-1")
	main_2 <<- annotate_go_terms("main-2")
	main_3 <<- annotate_go_terms("main-3")
	main_4 <<- annotate_go_terms("main-4")
	main_5 <<- annotate_go_terms("main-5")
	main_6 <<- annotate_go_terms("main-6")
	main_7 <<- annotate_go_terms("main-7")
	main_8 <<- annotate_go_terms("main-8")

	main_1_1 <<- annotate_go_terms("main-1-1")
	main_1_2 <<- annotate_go_terms("main-1-2")
	main_1_3 <<- annotate_go_terms("main-1-3")
	main_1_4 <<- annotate_go_terms("main-1-4")
	main_1_5 <<- annotate_go_terms("main-1-5")

	main_2_1 <<- annotate_go_terms("main-2-1")
	main_2_2 <<- annotate_go_terms("main-2-2")
	main_2_3 <<- annotate_go_terms("main-2-3")
	main_2_4 <<- annotate_go_terms("main-2-4")

	main_3_1 <<- annotate_go_terms("main-3-1")
	main_3_2 <<- annotate_go_terms("main-3-2")
	main_3_3 <<- annotate_go_terms("main-3-3")

	main_4_1 <<- annotate_go_terms("main-4-1")
	main_4_2 <<- annotate_go_terms("main-4-2")

	main_5_1 <<- annotate_go_terms("main-5-1")
	main_5_2 <<- annotate_go_terms("main-5-2")
	main_5_3 <<- annotate_go_terms("main-5-3")
	main_5_4 <<- annotate_go_terms("main-5-4")
	main_5_5 <<- annotate_go_terms("main-5-5")

	main_6_1 <<- annotate_go_terms("main-6-1")
	main_6_2 <<- annotate_go_terms("main-6-2")

	main_7_1 <<- annotate_go_terms("main-7-1")
	main_7_2 <<- annotate_go_terms("main-7-2")
	main_7_3 <<- annotate_go_terms("main-7-3")
	main_7_4 <<- annotate_go_terms("main-7-4")

	main_8_1 <<- annotate_go_terms("main-8-1")
	main_8_2 <<- annotate_go_terms("main-8-2")
	main_8_3 <<- annotate_go_terms("main-8-3")

	plot_go_terms(main_1, "main-1", 10)
	plot_go_terms(main_2, "main-2", 10)
	plot_go_terms(main_3, "main-3", 10)
	plot_go_terms(main_4, "main-4", 10)
	plot_go_terms(main_5, "main-5", 10)
	plot_go_terms(main_6, "main-6", 10)
	plot_go_terms(main_7, "main-7", 10)
	plot_go_terms(main_8, "main-8", 10)


	plot_go_terms(main_1_1, "main-1-1", 10)
	plot_go_terms(main_1_2, "main-1-2", 10)
	plot_go_terms(main_1_3, "main-1-3", 10)
	plot_go_terms(main_1_4, "main-1-4", 10)
	plot_go_terms(main_1_5, "main-1-5", 10)

	plot_go_terms(main_2_1, "main-2-1", 10)
	plot_go_terms(main_2_2, "main-2-2", 10)
	plot_go_terms(main_2_3, "main-2-3", 10)
	plot_go_terms(main_2_4, "main-2-4", 10)

	plot_go_terms(main_3_1, "main-3-1", 10)
	plot_go_terms(main_3_2, "main-3-2", 10)
	plot_go_terms(main_3_3, "main-3-3", 10)

	plot_go_terms(main_4_1, "main-4-1", 10)
	plot_go_terms(main_4_2, "main-4-2", 10)

	plot_go_terms(main_5_1, "main-5-1", 10)
	plot_go_terms(main_5_2, "main-5-2", 10)
	plot_go_terms(main_5_3, "main-5-3", 10)
	plot_go_terms(main_5_4, "main-5-4", 10)
	plot_go_terms(main_5_5, "main-5-5", 10)

	plot_go_terms(main_6_1, "main-6-1", 10)
	plot_go_terms(main_6_2, "main-6-2", 10)

	plot_go_terms(main_7_1, "main-7-1", 10)
	plot_go_terms(main_7_2, "main-7-2", 10)
	plot_go_terms(main_7_3, "main-7-3", 10)
	plot_go_terms(main_7_4, "main-7-4", 10)

	plot_go_terms(main_8_1, "main-8-1", 10)
	plot_go_terms(main_8_2, "main-8-2", 10)
	plot_go_terms(main_8_3, "main-8-3", 10)
}

plot_go_clust_wt <- function() {
	a66_1 <<- annotate_go_terms("a66-1")
	a66_2 <<- annotate_go_terms("a66-2")
	a66_3 <<- annotate_go_terms("a66-3")
	a66_4 <<- annotate_go_terms("a66-4")
	a66_5 <<- annotate_go_terms("a66-5")
	a66_6 <<- annotate_go_terms("a66-6")

	a66_1_1 <<- annotate_go_terms("a66-1-1")
	a66_1_2 <<- annotate_go_terms("a66-1-2")

	a66_2_1 <<- annotate_go_terms("a66-2-1")
	a66_2_2 <<- annotate_go_terms("a66-2-2")
	a66_2_3 <<- annotate_go_terms("a66-2-3")
	a66_2_4 <<- annotate_go_terms("a66-2-4")

	a66_3_1 <<- annotate_go_terms("a66-3-1")
	a66_3_2 <<- annotate_go_terms("a66-3-2")
	a66_3_3 <<- annotate_go_terms("a66-3-3")
	a66_3_4 <<- annotate_go_terms("a66-3-4")

	a66_4_1 <<- annotate_go_terms("a66-4-1")
	a66_4_2 <<- annotate_go_terms("a66-4-2")

	a66_5_1 <<- annotate_go_terms("a66-5-1")
	a66_5_2 <<- annotate_go_terms("a66-5-2")
	# a66_5_3 <<- annotate_go_terms("a66-5-3")
	# a66_5_4 <<- annotate_go_terms("a66-5-4")
	# a66_5_5 <<- annotate_go_terms("a66-5-5")

	a66_6_1 <<- annotate_go_terms("a66-6-1")
	a66_6_2 <<- annotate_go_terms("a66-6-2")

	plot_go_terms(a66_1, "a66-1", 10)
	plot_go_terms(a66_2, "a66-2", 10)
	plot_go_terms(a66_3, "a66-3", 10)
	plot_go_terms(a66_4, "a66-4", 10)
	plot_go_terms(a66_5, "a66-5", 10)
	plot_go_terms(a66_6, "a66-6", 10)

	plot_go_terms(a66_1_1, "a66-1-1", 10)
	plot_go_terms(a66_1_2, "a66-1-2", 10)
	plot_go_terms(a66_2_1, "a66-2-1", 10)
	plot_go_terms(a66_2_2, "a66-2-2", 10)
	plot_go_terms(a66_2_3, "a66-2-3", 10)
	plot_go_terms(a66_2_4, "a66-2-4", 10)
	plot_go_terms(a66_3_1, "a66-3-1", 10)
	plot_go_terms(a66_3_2, "a66-3-2", 10)
	plot_go_terms(a66_3_3, "a66-3-3", 10)
	plot_go_terms(a66_3_4, "a66-3-4", 10)
	plot_go_terms(a66_4_1, "a66-4-1", 10)
	plot_go_terms(a66_4_2, "a66-4-2", 10)
	plot_go_terms(a66_5_1, "a66-5-1", 10)
	plot_go_terms(a66_5_2, "a66-5-2", 10)
	# plot_go_terms(a66_5_3, "a66-5-3", 10)
	# plot_go_terms(a66_5_4, "a66-5-4", 10)
	# plot_go_terms(a66_5_5, "a66-5-5", 10)
	plot_go_terms(a66_6_1, "a66-6-1", 10)
	plot_go_terms(a66_6_2, "a66-6-2", 10)
}


quickgo_clust_mut <- function() {
	quickgo1(main_1, "main-1", 20)
	quickgo1(main_2, "main-2", 20)
	quickgo1(main_3, "main-3", 20)
	quickgo1(main_4, "main-4", 20)
	quickgo1(main_5, "main-5", 20)
	quickgo1(main_6, "main-6", 20)
	quickgo1(main_7, "main-7", 20)
	quickgo1(main_8, "main-8", 20)

	quickgo1(main_1_1, "main-1-1", 20)
	quickgo1(main_1_2, "main-1-2", 20)
	quickgo1(main_1_3, "main-1-3", 20)
	quickgo1(main_1_4, "main-1-4", 20)
	quickgo1(main_1_5, "main-1-5", 20)

	quickgo1(main_2_1, "main-2-1", 20)
	quickgo1(main_2_2, "main-2-2", 20)
	quickgo1(main_2_3, "main-2-3", 20)
	quickgo1(main_2_4, "main-2-4", 20)

	quickgo1(main_3_1, "main-3-1", 20)
	quickgo1(main_3_2, "main-3-2", 20)
	quickgo1(main_3_3, "main-3-3", 20)

	quickgo1(main_4_1, "main-4-1", 20)
	quickgo1(main_4_2, "main-4-2", 20)

	quickgo1(main_5_1, "main-5-1", 20)
	quickgo1(main_5_2, "main-5-2", 20)
	quickgo1(main_5_3, "main-5-3", 20)
	quickgo1(main_5_4, "main-5-4", 20)
	quickgo1(main_5_5, "main-5-5", 20)

	quickgo1(main_6_1, "main-6-1", 20)
	quickgo1(main_6_2, "main-6-2", 20)

	quickgo1(main_7_1, "main-7-1", 20)
	quickgo1(main_7_2, "main-7-2", 20)
	quickgo1(main_7_3, "main-7-3", 20)
	quickgo1(main_7_4, "main-7-4", 20)

	quickgo1(main_8_1, "main-8-1", 20)
	quickgo1(main_8_2, "main-8-2", 20)
	quickgo1(main_8_3, "main-8-3", 20)

	quickgo2(main_1, main_2, "main-1-2", 10)

	quickgo2(main_1_1, main_1_2, "main-1-1-1-2", 10)
	quickgo2(main_1_3, main_1_2, "main-1-3-1-2", 8)
	quickgo2(main_1_4, main_1_2, "main-1-4-1-2", 10)
	quickgo2(main_1_5, main_1_2, "main-1-5-1-2", 10)

	quickgo2(main_1_2, main_2_2, "main-1-2-2-2", 10)
	quickgo2(main_2_1, main_2_2, "main-2-1-2-2", 8)
	quickgo2(main_2_3, main_2_2, "main-2-3-2-2", 10)
	quickgo2(main_2_4, main_2_2, "main-2-4-2-2", 10)

	quickgo2(main_8_1, main_8_3, "main-8-1-8-3", 10)
	quickgo2(main_4_1, main_4_2, "main-4-1-4-2", 10)
}

quickgo_clust_wt <- function() {
	quickgo1(a66_1, "a66-1", 20)
	quickgo1(a66_2, "a66-2", 20)
	quickgo1(a66_3, "a66-3", 20)
	quickgo1(a66_4, "a66-4", 20)
	quickgo1(a66_5, "a66-5", 20)
	quickgo1(a66_6, "a66-6", 20)

	quickgo1(a66_1_1, "a66-1-1", 20)
	quickgo1(a66_1_2, "a66-1-2", 20)
	quickgo1(a66_2_1, "a66-2-1", 20)
	quickgo1(a66_2_2, "a66-2-2", 20)
	quickgo1(a66_2_3, "a66-2-3", 20)
	quickgo1(a66_2_4, "a66-2-4", 20)
	quickgo1(a66_3_1, "a66-3-1", 20)
	quickgo1(a66_3_2, "a66-3-2", 20)
	quickgo1(a66_3_3, "a66-3-3", 20)
	quickgo1(a66_3_4, "a66-3-4", 20)
	quickgo1(a66_4_1, "a66-4-1", 20)
	quickgo1(a66_4_2, "a66-4-2", 20)
	quickgo1(a66_5_1, "a66-5-1", 20)
	quickgo1(a66_5_2, "a66-5-2", 20)
	quickgo1(a66_6_1, "a66-6-1", 20)
	quickgo1(a66_6_2, "a66-6-2", 20)

	quickgo2(a66_1, a66_2, "a66-1-2", 20)

	quickgo2(a66_1_1, a66_1_2, "a66-1-1-1-2", 20)
	quickgo2(a66_2_1, a66_2_2, "a66-2-1-2-2", 10)

	quickgo2(a66_3_1, a66_3_2, "a66-3-1-3-2", 10)

	quickgo2(a66_5_1, a66_5_2, "a66-5-1-5-2", 10)


	quickgo2(a66_1_1, a66_1_2, "a66-1-1-1-2", 10)
	quickgo2(a66_1_3, a66_1_2, "a66-1-3-1-2", 8)
	quickgo2(a66_1_4, a66_1_2, "a66-1-4-1-2", 10)
	quickgo2(a66_1_5, a66_1_2, "a66-1-5-1-2", 10)

	quickgo2(a66_1_2, a66_2_2, "a66-1-2-2-2", 10)
	quickgo2(a66_2_1, a66_2_2, "a66-2-1-2-2", 8)
	quickgo2(a66_2_3, a66_2_2, "a66-2-3-2-2", 10)
	quickgo2(a66_2_4, a66_2_2, "a66-2-4-2-2", 10)
}

go_genes_for_vero <- function(name) {
	d <- read.csv(paste0("../pip3-rna-seq-output/GO/", name, "/REVIGO.csv"))
	# d <- d[d$eliminated == 0, ]
	d <- d[,c(1, 2, 7)]
	d <- d[order(d$log10.p.value),]
	d <- as.data.table(d)

	t <- readRDS(paste0("../pip3-rna-seq-output/GO/", name, "/genes-ann-BP.rds"))
	t <- t[names(t) %in% d$term_ID]
	t1 <- lapply(seq_along(t),
				 function(y, n, i) {cbind(y[[i]], d[term_ID == n[[i]]])},
				 y = t, n = names(t))
	t2 <- do.call("rbind", t1)
	colnames(t2)[1] <- "ensembl_gene_id"
	setkey(t2, "ensembl_gene_id")
	mart <- readRDS("../pip3-rna-seq-output/rds/count-matrix-ann.RDS")
	mart <- as.data.table(mart)
	setkey(mart, "ensembl_gene_id")
	t3 <- t2[mart, allow.cartesian=TRUE]
	t3 <- t3[!is.na(term_ID)]
	setkey(t3, "term_ID")
	t3 <- t3[order(log10.p.value)]
	write.csv(t3, file = paste0("../pip3-rna-seq-output/data-for-vero/GO/", name, ".csv"), row.names = F)
}

butterfly_paper_comparisons <- function() {
        # analysis of similarities of our results with the results of a new
        # paper: Hart, J. R. et al. The butterfly effect in cancer:
        # A single base mutation can remodel the cell.
        # Proc. Natl. Acad. Sci. U. S. A. 112, 1131–1136 (2015).
        d <- read.csv("../pip3-rna-seq-input/GSE63452_mcf10a.vs.pik3ca.h1047r.csv", sep = ",")
        # select only significant genes at 0hr time point: 3485 genes
        d1 <- d[d$X0hr.p.value < 0.05, ]
        d1 <- hgnc_symbol_to_ensembl_id(d1$gene)
        initialize_sets()
        venn(list("Butterfly" = d1$ensembl_gene_id, "Our KI" = ki, "Our PTEN" = pten),
             TRUE, "comparison-with-new-paper-0hr")
        # correlations between WTs and KIs
        count.matrix <- readRDS("../pip3-rna-seq-output/rds/count-matrix.rds")
        c <- count.matrix[,c(1:3, 58:60)]
        c <- as.data.frame(c)
        c$ensembl_gene_id <- rownames(c)
        d2 <- d[,1:7]
        gene.names <- hgnc_symbol_to_ensembl_id(d2$gene)
        colnames(gene.names)[2] <- "gene"
        d2 <- merge(d2, gene.names)
        d2 <- merge(d2, c)
        cor_heatmap(as.matrix(d2[,c(3:14)]), "test.pdf")
        # plot a correlation matrix from a count matrix
        # calculate pearson's correlation coefficients
        cor.matrix <- cor(as.matrix(d2[,c(3:14)]), method = "pearson")
        # plot correlation matrix in a file with 'name'
        pdf(file = "../pip3-rna-seq-output/figures/cor-butterfly.pdf", w = 6, h = 6)
        heatmap.2(cor.matrix, Rowv = FALSE, Colv = FALSE, dendrogram = "none",
                  col=bluered(99), breaks = 100, trace = "none", keysize = 1.5, margins = c(10, 10))
        dev.off()
        
        res <- prcomp(as.matrix(d2[,c(3:14)]))
        pdf(file = "../pip3-rna-seq-output/figures/pca-variances-butterfly.pdf", w=7, h=6)
        print(screeplot(res))
        dev.off()
        data <- as.data.frame(res$rotation[,1:3])
        
        data$Condition <- c("KI", "KI", "KI", "WT", "WT", "WT", "KI", "KI", "KI", "WT", "WT", "WT")
        data$Paper <- c(rep("Butterfly", 6), rep("Ours", 6))
        p <- ggplot(data, aes(PC1,PC2, color = Condition)) +
                geom_point(aes(shape = Paper)) +
                theme_bw()
        pdf(file = "../pip3-rna-seq-output/figures/pca12-butterfly.pdf", w=5, h=4)
        print(p)
        dev.off()
        p <- ggplot(data, aes(PC1,PC3, color = Condition)) +
                geom_point(aes(shape = Paper)) +
                theme_bw()
        pdf(file = "../pip3-rna-seq-output/figures/pca13-butterfly.pdf", w=5, h=4)
        print(p)
        dev.off()
        p <- ggplot(data, aes(PC2,PC3, color = Condition)) +
                geom_point(aes(shape = Paper)) +
                theme_bw()
        pdf(file = "../pip3-rna-seq-output/figures/pca23-butterfly.pdf", w=5, h=4)
        print(p)
        dev.off()
}

prepare_data_for_len_phil_effect <- function(){
        venn(list("PTEN" = pten, "KI" = ki, "A66" = a66.tc.cond), T, "len-phil-effect")
        l.p.genes <- pten[pten %in% ki & pten %in% a66.tc.cond]
        clust_boot(l.p.genes, 2, 8, "len-phil-effect", 0)
        plot_bootstrap_data("len-phil-effect")
        plot_genes(l.p.genes, F, "len-phil-effect")
        clust <- get_clust_genes("len-phil-effect", "0-2")
        clusts <- list(clust$partition)
        plot_all_clusts(clusts, "len-phil-effect", c(0))
        ind <- 1
        for(j in clusts) {
                for (i in 1:length(unique(j))) {
                        GO(names(j[j == i]), all.genes, paste0("len-phil-effect", c(1, 2)[ind], "-", i), 0.05)
                }
                ind <- ind + 1
        }
        GO(l.p.genes, all.genes, "len-phil-effect", 0.05)
        plot_genes(names(clusts[[1]][clusts[[1]] == 1]), F, "len-phil-effect-1")
        plot_genes(names(clusts[[1]][clusts[[1]] == 2]), F, "len-phil-effect-2")
        
        len_phil_effect1 <- annotate_go_terms("len-phil-effect1-1")
        len_phil_effect2 <- annotate_go_terms("len-phil-effect1-2")
        len_phil_effect <- annotate_go_terms("len-phil-effect")
        plot_go_terms(len_phil_effect1, "len-phil-effect1", 10)
        plot_go_terms(len_phil_effect2, "len-phil-effect2", 10)
        plot_go_terms(len_phil_effect, "len-phil-effect", 10)
        
        d <- read.csv("../pip3-rna-seq-input/GSE63452_mcf10a.vs.pik3ca.h1047r.csv", sep = ",")
        # select only significant genes at 0hr time point: 3485 genes
        d1 <- d[d$X0hr.p.value < 0.05, ]
        d1 <- hgnc_symbol_to_ensembl_id(d1$gene)
        venn(list("Butterfly" = d1$ensembl_gene_id, "'Expected' Genes" = l.p.genes),
             TRUE, "comparison-with-new-paper-expected-effect")
        plot_genes(l.p.genes[l.p.genes %in% d1$ensembl_gene_id], F, "len-phil-effect-in-butterfly")
        plot_genes(l.p.genes[!l.p.genes %in% d1$ensembl_gene_id], F, "len-phil-effect-not-in-butterfly")
}
