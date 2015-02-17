
linc_rna <- function(cond) {
	gene.ann <- read.csv("../pip3-rna-seq-input/linc-rna/from-biomart-11-07-14.txt")
	if(cond == "mut") {
	        res <- gene.ann[gene.ann[,1] %in% ki[ki %in% pten], ]
	} else {
	        res <- gene.ann[gene.ann[,1] %in% a66_nost_eff, ] 
	}

	lincRNA <- res[res[,3] == "lincRNA", ]
	genes <- unique(lincRNA[,1])
	plot_genes(genes, F, paste0("lincRNA-", cond))
	plot_genes(genes, T, paste0("lincRNA-", cond))
}

plot_linc_rna <- function() {
	files <- list.files("../pip3-rna-seq-input/linc-rna/", full.name = T)
	files <- files[grepl(".txt", files)]
	files <- files[!grepl("biomart", files)]

	for(file in files) {
		d <- read.table(file, header = F)
		name <- strsplit(file, "\\/")[[1]][5]
		name <- strsplit(name, "\\.")[[1]][1]
		plot_genes(d[,1], F, paste0("linc-rna-", name, "-raw"))
		plot_genes(d[,1], T, paste0("linc-rna-", name, "-norm"))
	}
}
