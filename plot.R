ensembl_id_to_hgnc_symbol <- function(genes) {
        mart <- readRDS("../pip3-rna-seq-output/rds/count-matrix-ann.RDS")
        res <- unique(mart[mart$ensembl_gene_id %in% genes, 1:2])
        res <- res[res$hgnc_symbol != "", ]
        return(res)
}

hgnc_symbol_to_ensembl_id <- function(genes) {
        mart <- readRDS("../pip3-rna-seq-output/rds/count-matrix-ann.RDS")
        res <- unique(mart[mart$hgnc_symbol %in% genes, 1:2])
        # res <- res[res$ensembl_gene_id != "", ]
        return(res)
}

get_plot_data <- function(genes, norm) {
        # if genes are represented by ensembl ids
        if(grepl("ENSG000", genes[1])){
        res <- ensembl_id_to_hgnc_symbol(genes)
        } else {
        res <- hgnc_symbol_to_ensembl_id(genes)
        }
        colnames(res) <- c("id", "name")
        
        if(norm) {
        plot.data <- readRDS("../pip3-rna-seq-output/rds/plot-time-courses-all-genes-norm-by-length.rds")
        } else {
        plot.data <- readRDS("../pip3-rna-seq-output/rds/plot-time-courses-all-genes.rds")
        }
        plot.data <- plot.data[plot.data$id %in% res$id, ]
        plot.data <- merge(plot.data, res)
        return(list(plot.data, res))
}

plot_genes <- function(genes, norm, s) {
        plot.data <- get_plot_data(genes, norm)[[1]]
        res <- get_plot_data(genes, norm)[[2]]
        
        limits <- aes(ymax = ymax, ymin = ymin)
        
        p <- ggplot(plot.data, aes(time, value, group = cond, color = cond)) +
        geom_line(size = 1) +
        geom_point(size = 3) +
        facet_wrap(~ name, scale = "free_y") +
        geom_errorbar(limits, size = 0.5, width = 5) +
        labs(x = "Time, min", y = "Read counts") +
        theme_bw()
        if(norm) {
        file.name <- paste0("../pip3-rna-seq-output/figures/genes-norm-", s, ".pdf")
        } else {
        file.name <- paste0("../pip3-rna-seq-output/figures/genes-", s, ".pdf")
        }
        pdf(file = file.name,
        width = (sqrt(length(res$name)) + 2)*3,
        height = sqrt(length(res$name))*3)
        print(p)
        dev.off()
}

plot_prdm1_genes <- function(genes, norm, s) {
        plot.data <- get_plot_data(genes, norm)[[1]]
        res <- get_plot_data(genes, norm)[[2]]
        
        limits <- aes(ymax = ymax, ymin = ymin)
        
        plot.data$name <- factor(plot.data$name, levels = genes)
        
        p <- ggplot(plot.data,
        aes(time, value, group = cond, color = cond)) +
        geom_point(data = plot.data[plot.data$name %in% genes & plot.data$time == 300 & plot.data$cond == "a66_nost",]) +
        geom_line() + facet_grid(name~., scale = "free") +
        geom_errorbar(limits, width = 0.25) +
        theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text=element_text(size = 7),
            legend.position="none",
            plot.background=element_blank(),
            panel.background=element_blank(),
            strip.text = element_blank(),
            strip.background = element_blank())
        # theme_tufte()
        if(norm) {
        file.name <- paste0("../pip3-rna-seq-output/figures/genes-norm-", s, ".pdf")
        } else {
        file.name <- paste0("../pip3-rna-seq-output/figures/genes-", s, ".pdf")
        }
        pdf(file = file.name,
        width = 2,
        height = 10)
        print(p)
        dev.off()
}

# plot_prdm1_genes_multiple <- function(genes, norm, s) {
#   plot.data <- get_plot_data(genes, norm)[[1]]
#   res <- get_plot_data(genes, norm)[[2]]

#   limits <- aes(ymax = ymax, ymin = ymin)

#   for(g in genes) {
#     p <- ggplot(plot.data[plot.data$name == g,],
#       aes(time, value, group = cond, color = cond)) +
#         geom_point(data = plot.data[plot.data$name == g & plot.data$time == 300 & plot.data$cond == "a66_nost",], size = 3) +
#         geom_line(size = 1.5)+
#         geom_errorbar(limits, width = 0.25) +
#         theme(axis.title.x=element_blank(),
#               axis.title.y=element_blank(),
#               legend.position="none",
#               plot.background=element_blank(),
#               panel.background=element_blank())
#         # theme_bw()
#         # theme_tufte()
#     if(norm) {
#       file.name <- paste0("../pip3-rna-seq-output/figures/genes-PRDM1-norm-", g, ".png")
#     } else {
#       file.name <- paste0("../pip3-rna-seq-output/figures/genes-PRDM1-", g, ".png")
#     }
#     ggsave(p, file = file.name, dpi = 90)
#   }
# }

plot_genes_by_motif <- function(mot.table, norm, s) {
  mot.table <- mot.table[order(score, decreasing = T)]
  plot.data <- get_plot_data(mot.table[,target], norm)[[1]]
  mot.table$name <- mot.table$target
  plot.data <- merge(plot.data, mot.table[,list(motif, score, name)])

  limits <- aes(ymax = ymax, ymin = ymin)

  p <- ggplot(plot.data,
    aes(time, value, group = cond, color = cond)) +
      geom_line() +
      facet_grid(name ~ motif, scale = "free_y") +
      geom_errorbar(limits, width = 0.25)
  if(norm) {
    file.name <- paste0("../pip3-rna-seq-output/figures/genes-norm-", s, ".pdf")
  } else {
    file.name <- paste0("../pip3-rna-seq-output/figures/genes-", s, ".pdf")
  }
  pdf(file = file.name, width = 20, height = 20)
  print(p)
  dev.off()
}

plot_mirna <- function(mirnas, s) {
  # if genes are represented by ensembl ids
  # if(grepl("ENSG000", genes[1])){
  #   res <- ensembl_id_to_hgnc_symbol(genes)
  # } else {
  #   res <- hgnc_symbol_to_ensembl_id(genes)
  # }
  # colnames(res) <- c("id", "name")
  
  plot.data <- readRDS("../pip3-rna-seq-output/rds/plot-time-courses-all-mirna.rds")
  plot.data <- plot.data[plot.data$id %in% mirnas, ]

  limits <- aes(ymax = ymax, ymin = ymin)

  p <- ggplot(plot.data,
    aes(time, value, group = cond, color = cond)) +
      geom_line() + facet_wrap(~ id, scale = "free_y") +
      geom_errorbar(limits, width = 0.25)
  pdf(file = paste0("../pip3-rna-seq-output/figures/mirna-", s, ".pdf"),
    width = (sqrt(length(mirnas)) + 2)*3,
    height = sqrt(length(mirnas))*3)
  print(p)
  dev.off()
}

plot_exons <- function(cond, genes) {
  d <- readRDS(paste0("../pip3-rna-seq-output/rds/exon_diff_expr_wt_0_", cond, "_0.rds"))

  # filter by genes
  groups <- unique(d[,1])
  groups <- groups[grepl(paste(genes, collapse = "|"), groups)]
  d <- d[d[,1] %in% groups, ]

  # filter by padj
  groups <- unique(d[!is.na(d$padj) & d$padj < 0.01, 1])
  d <- d[d[,1] %in% groups, ]

  unlink(paste0("../rna-seq-media/DEXSeqReport_", cond), recursive = T)
  DEXSeqHTML( d, path = paste0("../rna-seq-media/DEXSeqReport_", cond), FDR = 0.01,
    color=c("#FF000080", "#0000FF80") )
}

plot_genes_by_constitutive_eff <- function(genes, s) {
  # query biomart
  if(grepl("ENSG000", genes[1, 1])){
    res <- ensembl_id_to_hgnc_symbol(genes[, 1])
  } else {
    res <- hgnc_symbol_to_ensembl_id(genes[, 1])
  }
  colnames(res) <- c("id", "name")

  res <- merge(res, genes)

  plot.data <- readRDS("../pip3-rna-seq-output/rds/plot-time-courses-all-genes.rds")
  plot.data <- plot.data[plot.data$id %in% res$id, ]
  plot.data <- merge(plot.data, res)

  plot.data <- plot.data[order(-plot.data$log2F.eff),]

  plot.data$name <- factor(plot.data$name, levels = unique(plot.data$name))

  limits <- aes(ymax = ymax, ymin = ymin)

  p <- ggplot(plot.data,
    aes(time, value, group = cond, color = cond)) +
      geom_line() + facet_wrap(~ name, scale = "free_y") +
      geom_errorbar(limits, width = 0.25)
  pdf(file = paste0("../pip3-rna-seq-output/figures/genes-", s, ".pdf"),
    width = (sqrt(length(res$name)) + 2)*3,
    height = sqrt(length(res$name))*3)
  print(p)
  dev.off()
}

plot_interesting_genes <- function() {
  plot_genes(c("ACPL2", "ATP9A", "AXL", "BCL7A", "CLASP1", "GLG1", "GPC1",
    "KRT10", "LPAR3", "NTN4", "RGCC", "SDF2L1", "THBD"),
    "A66+KI+PTEN+EGF")

  plot_genes(c("CCAT1", "DAW1", "PI3", "ABCC2", "TACC2", "TRIML2", "NRG4", "TNNI2",
    "ZNF584", "IL27RA", "EPGN", "STK40", "PIK3R1", "PROM2", "BAK1", "VAV2",
    "CAPN6", "VAV1", "KRT13", "MYEOV", "PLA2G4F", "CES1", "AKR1B15", "PTK6",
    "S100A8", "PREX1", "H19", "S100A9", "SYT8", "NQO1", "HMGN5", "MED24", "ATP5G3"),
    "KI+PTEN+EGF_NO A66_UP-UP")

  plot_genes(c("NGFRAP1", "DTX3", "LBH", "ATP8B2", "ID4", "COL4A3", "CD27-AS1",
    "CALHM2", "MYL9", "ASPA", "VWA5A", "PXDN", "TMEM139", "STON1", "FRMD3",
    "LMCD1", "CLDN1", "COL4A4", "LPXN", "AKAP12", "F2R", "PIK3IP1", "ARHGAP29",
    "CAMK2N1", "INPP4B", "LIF", "SLC46A3", "KIF3C", "DENND2D", "NNMT",
    "ARHGEF28", "REPS2", "HPSE", "SLFN5", "HEG1", "CREB3L4", "C16orf62",
    "TCP11L1", "PTK2", "RASSF5", "RNF152", "TPM1", "HIP1", "CALCOCO1",
    "SWAP70", "CAMLG", "SPTBN1", "PDLIM1", "HERPUD2"),
    "KI+PTEN+EGF_NO A66_ DOWN-DOWN")

  plot_genes(c("RGS11", "SPARC", "CDH13", "SRPX", "C8orf47", "QPCT", "NME5",
    "MAGEH1", "SLIT3", "ZNF37", "BDH1", "THNSL2", "MCAM", "PLA2R1", "ZNF585A",
    "CFI", "ZNF358", "FN1", "ITGB2", "NMNAT2", "NMNAT3", "DOCK10", "DAB2",
    "ATL1", "TUBA1A", "LOX", "ARG2", "CARD16", "SCEL", "FAP", "FSTL1",
    "LTBP1", "PKN1", "DIXDC1", "CRIM1", "AKT3", "APLP2"),
    "KI+PTEN_NO A66_ NO EGF _DOWN-DOWN")
}

group_genes_by_mean_expr <- function(genes) {
  count.matrix.av <- readRDS("../rna-seq-media/data/count_matrix_av.RDS")
  rownames(count.matrix.av) <- count.matrix.av$id
  wt <- as.matrix(count.matrix.av[, 2:7])
  a66 <- as.matrix(count.matrix.av[, 8:13])
  pten <- as.matrix(count.matrix.av[, 14:19])
  ki <- as.matrix(count.matrix.av[, 20:25])
  correlated.pten <- (pten + 1)/(wt + 1)/6
  correlated.a66 <- (a66 + 1)/(wt + 1)/6
  correlated.ki <- (ki + 1)/(wt + 1)/6
  count.matrix.av$cor.pten <- log(rowSums(correlated.pten))
  count.matrix.av$cor.a66 <- log(rowSums(correlated.a66))
  count.matrix.av$cor.ki <- log(rowSums(correlated.ki))

  test <- data.frame(log2.A66.WT = count.matrix.av$cor.a66,
    log2.PTEN.WT = count.matrix.av$cor.pten, log2.ki.WT = count.matrix.av$cor.ki)
  test <- melt(test)
  p <- ggplot(test, aes(value, fill = variable)) +
    geom_density(aes(color = variable), alpha = 0.2) +
    coord_cartesian(xlim = c(-2, 2))
  pdf(file = "../rna-seq-media/figures/expr_ratios.pdf")
  print(p)
  dev.off()

  test <- data.table(data.frame(id = count.matrix.av$id, a66 = count.matrix.av$cor.a66,
    pten = count.matrix.av$cor.pten, ki = count.matrix.av$cor.ki))

  test <- test[id %in% genes]

  eff1 <- as.character(test[pten > 0.2 & ki > 0.2 & a66 < -0.1, id])
  eff2 <- as.character(test[pten > 0.2 & ki > 0.2 & a66 > 0.1, id])
  eff3 <- as.character(test[pten > 0.2 & ki < -0.2 & a66 < -0.1, id])
  eff4 <- as.character(test[pten > 0.2 & ki < -0.2 & a66 > 0.1, id])
  eff5 <- as.character(test[pten < -0.2 & ki > 0.2 & a66 < -0.1, id])
  eff6 <- as.character(test[pten < -0.2 & ki > 0.2 & a66 > 0.1, id])
  eff7 <- as.character(test[pten < -0.2 & ki < -0.2 & a66 < -0.1, id])
  eff8 <- as.character(test[pten < -0.2 & ki < -0.2 & a66 > 0.1, id])
  if(length(eff1 > 0)) {
    plot_genes(eff1, "eff1")
  }
  if(length(eff2 > 0)) {
    plot_genes(eff2, "eff2")
  }
  if(length(eff3 > 0)) {
    plot_genes(eff3, "eff3")
  }
  if(length(eff4 > 0)) {
    plot_genes(eff4, "eff4")
  }
  if(length(eff5 > 0)) {
    plot_genes(eff5, "eff5")
  }
  if(length(eff6 > 0)) {
    plot_genes(eff6, "eff6")
  }
  if(length(eff7 > 0)) {
    plot_genes(eff7, "eff7")
  }
  if(length(eff8 > 0)) {
    plot_genes(eff8, "eff8")
  }

  return(list(eff1 = eff1, eff2 = eff2, eff3 = eff3, eff4 = eff4, eff5 = eff5,
    eff6 = eff6, eff7 = eff7, eff8 = eff8))
}
