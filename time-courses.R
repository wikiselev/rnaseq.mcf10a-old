######################
# time course analysis
######################

order_time_courses_by_wt_peaks <- function() {
  c <- readRDS("../pip3-rna-seq-output/rds/count-matrix-av-scaled.rds")
  rownames(c) <- c$id
  c <- as.matrix(c[ , c(2:26)])
  maxs <- rowMax(c[ , c(1:6)])
  t <- as.data.frame(c)
  t$max <- maxs
  t$peak <- apply(t, 1, function(x) match(x[26], x[1:6]))

  egf.wt.resp <- c(egf.wt15, egf.wt40, egf.wt90, egf.wt180, egf.wt300)

  expr <- readRDS("../pip3-rna-seq-output/rds/count-matrix-av.rds")
  t$change15 <- expr[,3]/expr[,2]
  t$change40 <- expr[,4]/expr[,2]
  t$change90 <- expr[,5]/expr[,2]
  t$change180 <- expr[,6]/expr[,2]
  t$change300 <- expr[,7]/expr[,2]

  g2 <- rownames(t[rownames(t) %in% egf.wt.resp & t$peak == 2 & t$change15 > 1.5, ])
  g3 <- rownames(t[rownames(t) %in% egf.wt.resp & t$peak == 3 & t$change40 > 1.5, ])
  g4 <- rownames(t[rownames(t) %in% egf.wt.resp & t$peak == 4 & t$change90 > 2, ])
  g5 <- rownames(t[rownames(t) %in% egf.wt.resp & t$peak == 5 & t$change180 > 2, ])
  g6 <- rownames(t[rownames(t) %in% egf.wt.resp & t$peak == 6 & t$change300 > 2, ])

  t <- t[rownames(t) %in% c(g2, g3, g4, g5, g6), ]
  t$id <- rownames(t)

  names <- ensembl_id_to_hgnc_symbol(rownames(t))
  colnames(names) <- c("id", "gene_name")

  t <- merge(t, names)
  rownames(t) <- t$id  
  t <- t[,2:34]

  t1 <- t
  rownames(t1) <- t1$gene_name
  t1 <- t1[,1:32]
  t1 <- t1[order(t1$peak),]
  t1 <- as.matrix(t1)

  a66.names1 <- ensembl_id_to_hgnc_symbol(a66.set1)

  pdf(file = paste0("../pip3-rna-seq-output/figures/peaks-wt-23-!a661.pdf"),
    w = 10, h = 10)
  heatmap.2(t1[!(rownames(t1) %in% a66.names1[,2]) & (t1[,27] == 2 | t1[,27] == 3), c(1:24)], Rowv = TRUE, Colv = FALSE, dendrogram = "row",
    col=bluered(19), breaks = 20, trace = "none", keysize = 0.5)
  dev.off()

  pdf(file = paste0("../pip3-rna-seq-output/figures/peaks-wt-23-a661.pdf"),
    w = 10, h = 10)
  heatmap.2(t1[(rownames(t1) %in% a66.names1[,2]) & (t1[,27] == 2 | t1[,27] == 3), c(1:24)], Rowv = TRUE, Colv = FALSE, dendrogram = "row",
    col=bluered(19), breaks = 20, trace = "none", keysize = 0.5)
  dev.off()

  pdf(file = paste0("../pip3-rna-seq-output/figures/peaks-wt-4-!a661.pdf"),
    w = 10, h = 16)
  heatmap.2(t1[!(rownames(t1) %in% a66.names1[,2]) & t1[,27] == 4, c(1:24)], Rowv = TRUE, Colv = FALSE, dendrogram = "row",
    col=bluered(19), breaks = 20, trace = "none", keysize = 0.5)
  dev.off()

  pdf(file = paste0("../pip3-rna-seq-output/figures/peaks-wt-4-a66.pdf"),
    w = 10, h = 10)
  heatmap.2(t1[(rownames(t1) %in% a66.names1[,2]) & t1[,27] == 4, c(1:24)], Rowv = TRUE, Colv = FALSE, dendrogram = "row",
    col=bluered(19), breaks = 20, trace = "none", keysize = 0.5)
  dev.off()

  # pdf(file = paste0("../pip3-rna-seq-output/figures/peaks-wt-5.pdf"),
  #   w = 10, h = 45)
  # heatmap.2(t1[t1[,27] == 5, c(1:24)], Rowv = TRUE, Colv = FALSE, dendrogram = "row",
  #   col=bluered(19), breaks = 20, trace = "none", keysize = 0.5)
  # dev.off()

  # pdf(file = paste0("../pip3-rna-seq-output/figures/peaks-wt-6.pdf"),
  #   w = 10, h = 35)
  # heatmap.2(t1[t1[,27] == 6, c(1:24)], Rowv = TRUE, Colv = FALSE, dendrogram = "row",
  #   col=bluered(19), breaks = 20, trace = "none", keysize = 0.5)
  # dev.off()

  saveRDS(t, "../pip3-rna-seq-output/rds/peaks-wt-ens.rds")
  saveRDS(t1, "../pip3-rna-seq-output/rds/peaks-wt-names.rds")  
}

plot_time_courses_by_wt_peaks <- function(genes, name) {
  t <- readRDS("../pip3-rna-seq-output/rds/peaks-wt.rds")
  t <- t[rownames(t) %in% genes, ]
  # wt
  pdf(file = paste0("../pip3-rna-seq-output/figures/peaks-wt.pdf"),
    w = 8, h = 11)
  heatmap.2(t[t[,26] == 3,1:6], Colv = F)
  dev.off()
  # a66
  pdf(file = paste0("../rna-seq-media/figures/a66_tc_by_wt_peaks_", name, ".pdf"),
    w = 8, h = 11)
  heatmap.3(t[,7:12], Colv = F, scale = c("row"))
  dev.off()
  # pten
  pdf(file = paste0("../rna-seq-media/figures/pten_tc_by_wt_peaks_", name, ".pdf"),
    w = 8, h = 11)
  heatmap.3(t[,13:18], Colv = F, scale = c("row"))
  dev.off()
  # ki
  pdf(file = paste0("../rna-seq-media/figures/ki_tc_by_wt_peaks_", name, ".pdf"),
    w = 8, h = 11)
  heatmap.3(t[,19:24], Colv = F, scale = c("row"))
  dev.off()
}

plot_peaks_old <- function() {
  d <- get_node_genes()
  genes <- d$ki[d$ki %in% d$pten]
  egf_genes <- get_wt_egf_genes()
  g <- d$ki[d$ki %in% d$pten]
  g <- g[g %in% egf_genes]
  plot_time_courses_by_wt_peaks(g, "ki&pten&egf_old")

  g <- d$ki[d$ki %in% d$pten]
  g <- g[!(g %in% egf_genes)]
  plot_time_courses_by_wt_peaks(g, "ki&pten_old")
}

plot_amplitudes_old <- function() {
  d <- get_node_genes()
  genes <- d$ki[d$ki %in% d$pten]
  egf_genes <- get_wt_egf_genes()
  g <- d$ki[d$ki %in% d$pten]
  g <- g[g %in% egf_genes]


  c <- readRDS("../rna-seq-media/data/count_matrix_av.RDS")
  rownames(c) <- c$id
  c <- as.matrix(c[ , c(2:7, 14:26)])

  heatmap.3(log(c[rownames(c) %in% g, ]+1), Rowv = TRUE, Colv = FALSE,
    cluster.by.row = TRUE, dendrogram = "row")

  g <- d$ki[d$ki %in% d$pten]
  g <- g[!(g %in% egf_genes)]

  heatmap.3(log(c[rownames(c) %in% g, ]+1), Rowv = TRUE, Colv = FALSE,
    cluster.by.row = TRUE, dendrogram = "row")
}

plot_differential_amplitude <- function() {
  d <- get_node_genes()
  genes <- d$ki[d$ki %in% d$pten]
  egf_genes <- get_wt_egf_genes()
  g <- d$ki[d$ki %in% d$pten]
  g <- g[g %in% egf_genes]


  c <- readRDS("../rna-seq-media/data/count_matrix_av.RDS")
  rownames(c) <- c$id

  wt <- as.matrix(c[ , c(2:7)])
  pten <- as.matrix(c[ , c(14:19)])
  ki <- as.matrix(c[ , c(20:25)])

  pten.dif <- (pten+1)/(wt+1)
  ki.dif <- (ki+1)/(wt+1)

  m <- cbind(pten.dif, ki.dif)
  m <- log(m[rownames(m) %in% g, ])

  m.min  <- min(m[ m!=0 ], na.rm=TRUE)
  m.max <- max( m, na.rm = TRUE )
  pairs.breaks <- c(0, seq( m.min, m.max, length.out=30) )

  breaks=seq(-m.max, m.max, by=0.1)
  breaks=append(breaks, m.min)
  # mycol <- colorpanel(n=length(breaks)-1,low="green",mid="black",high="red")


  heatmap.3(m, Rowv = TRUE, Colv = FALSE,
    cluster.by.row = TRUE, dendrogram = "row", color.FUN = "bluered",
    breaks=breaks)

  g <- d$ki[d$ki %in% d$pten]
  g <- g[!(g %in% egf_genes)]

  heatmap.3(log(m[rownames(m) %in% g, ]), Rowv = TRUE, Colv = FALSE,
    cluster.by.row = TRUE, dendrogram = "row")

}



plot_a66 <- function() {
  c <- readRDS("../rna-seq-media/data/count_matrix_av.RDS")
  rownames(c) <- c$id
  wt <- as.matrix(c[ , c(2:7)])
  a66 <- as.matrix(c[ , c(8:13)])
  a66.dif <- a66-wt

  d <- get_diff_expr("wt_ko", 0.01)
  a66 <- rownames(d$ko_vs_wt)
  t <- a66.dif[rownames(a66.dif) %in% a66, ]

  heatmap.3(t, Rowv = TRUE,
    Colv = FALSE, cluster.by.row = TRUE, dendrogram = "row", scale = c("row"))
}

######################
# time course analysis
######################
