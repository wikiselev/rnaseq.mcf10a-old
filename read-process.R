count_matrix_from_files <- function(files, normalized) {
  # this function create a count matrix from raw htseq read counts files
  # it uses DESeq2 library to construct the matrix

  # create a design matrix for DESeq2
  sample.name <- sapply(strsplit(files, "_trimmed"), "[[", 1)
  condition <- sapply(strsplit(sample.name, "_"), "[[", 1)
  samples <- data.frame(sample.name = sample.name,
    file.name = paste0("../pip3-rna-seq-input/htseq-read-count-raw/", files),
    condition = condition)
  samples$condition <- factor(samples$condition,
    levels = unique(samples$condition))

  # import data from files using 'samples' design matrix
  cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples,
    directory = ".", design = ~ condition)

  # if normalized is TRUE then the count matrix is normalized by library size
  # using estimateSizeFactors function of DESeq2
  if (normalized) {
    cds <- estimateSizeFactors(cds)
    return(counts(cds, normalized = TRUE))
  } # otherwise return the raw matrix
  else {return(counts(cds))}
}

count_matrix_norm_matrix <- function(matrix) {
  colData <- data.frame(condition = factor(colnames(matrix), levels = colnames(matrix)))
  dds <- DESeqDataSetFromMatrix(countData = matrix, colData = colData,
    design = ~ condition)
  dds <- estimateSizeFactors(dds)
  return(counts(dds, normalized = TRUE))
}

cor_heatmap <- function(count.matrix, name) {
  # plot a correlation matrix from a count matrix
  # calculate pearson's correlation coefficients
  cor.matrix <- cor(count.matrix, method = "pearson")

  # plot correlation matrix in a file with 'name'
  pdf(file = paste0("../pip3-rna-seq-output/figures/", name), w = 10, h = 10)
  heatmap.2(cor.matrix, Rowv = FALSE, Colv = FALSE, dendrogram = "none",
    col=bluered(9), breaks = 10, trace = "none")
  dev.off()
}

arrange_count_matrix_columns <- function(count.matrix) {
  # rearrage columns in count matrix so that they ordered by treatments and times
  # if count matrix is without bad replicates - there are 75 samples
  if (dim(count.matrix)[2] == 75) {
    x <- c(1,7,13,2,8,14,5,11,17,6,12,18,3,9,15,4,10,16,19,25,31,20,26,32,23,29,35,
      24,30,36,21,27,33,22,28,34,37,38,39,40,46,52,41,47,53,44,50,56,45,51,57,42,48,54,43,
      49,55,58,73,68,59,64,69,62,74,75,63,67,72,60,65,70,61,66,71)
  } else { # with bad replicates - there are 78 samples
    x <- c(1,7,13,2,8,14,5,11,17,6,12,18,3,9,15,4,10,16,19,25,31,20,26,32,23,29,35,
      24,30,36,21,27,33,22,28,34,37,38,39,40,46,52,41,47,53,44,50,56,45,51,57,42,48,54,43,
      49,55,58,64,70,59,65,71,62,68,74,63,69,75,60,66,72,61,67,73,76,77,78)
  }
  return(count.matrix[,x])
}

count_matrix_2_data_table <- function(m) {
  # this function create a data table from a read count matrix
  # data tables are useful for some analysis
  d <- as.data.frame(m)
  d$gene_id <- rownames(d)
  d <- melt(d)
  d <- as.data.table(d)
  d[,cond:=sapply(strsplit(as.character(d[,variable]), "_"), "[[", 1)]
  d[,rep:=sapply(strsplit(as.character(d[,variable]), "_"), "[[", 2)]
  d[,time:=sapply(strsplit(as.character(d[,variable]), "_"), "[[", 3)]
  d <- d[,list(gene_id,value,cond,rep,time)]
  # change wrong treatment names to proper ones
  d[cond == "ko", cond:="a66"]
  d[cond == "konost", cond:="a66_nost"]
  return(d)
}

plot_read_density <- function(d, name) {
  d$facet <- paste(d$cond, d$rep, d$time, sep = "_")
  p <- ggplot(as.data.frame(d), aes(value, color = facet)) + geom_density() +
    scale_x_log10()
    # facet_grid(cond ~ time)
  pdf(file = paste0("../pip3-rna-seq-output/figures/read-density-", name, ".pdf"),
    w = 12, h = 12)
  print(p)
  dev.off()
}

plot_read_ecdf <- function(d, name) {
  d$facet <- paste(d$cond, d$rep, d$time, sep = "_")
  p <- ggplot(as.data.frame(d), aes(value, color = facet)) + stat_ecdf() +
    scale_x_log10()
    # facet_grid(cond ~ time)
  pdf(file = paste0("../pip3-rna-seq-output/figures/read-ecdf-", name, ".pdf"),
    w = 12, h = 12)
  print(p)
  dev.off()
}

average_repl <- function(m, cond, time) {
  # average count matrix m over replicates at condition 'cond' and time 'time'
  filt <- sapply(paste0(cond, "_[0123456789]_", time), grepl, colnames(m))[,1]
  return(rowMeans(m[, filt]))
}

sd_repl <- function(m, cond, time) {
  # average count matrix m over replicates at condition 'cond' and time 'time'
  filt <- sapply(paste0(cond, "_[0123456789]_", time), grepl, colnames(m))[,1]
  return(apply(m[, filt], 1, sd))
}

melt_data <- function(m) {
  # melt data in create_plot_table function
  d <- melt(m)
  d <- as.data.table(d)
  d[,cond:=sapply(strsplit(as.character(d[,variable]), "_"), "[[", 1)]
  d[,time:=sapply(strsplit(as.character(d[,variable]), "_"), "[[", 2)]
  d <- d[,list(id,value,cond,time)]
  d[cond == "ko", cond:="a66"]
  d[cond == "konost", cond:="a66_nost"]
  d <- as.data.frame(d)
  d$time <- as.numeric(d$time)
  return(d)
}

average_count_matrix <- function(count.matrix, conditions) {
  times <- c(0, 15, 40, 90, 180, 300)
  # do the averaging
  count.matrix.av <- data.frame(id = rownames(count.matrix))
  count.matrix.sd <- data.frame(id = rownames(count.matrix))
  for(cond in conditions) {
    for(time in times) {
      count.matrix.av$new <- average_repl(count.matrix, cond, time)
      colnames(count.matrix.av)[dim(count.matrix.av)[2]] <- paste0(cond, "_", time)
      count.matrix.sd$new <- sd_repl(count.matrix, cond, time)
      colnames(count.matrix.sd)[dim(count.matrix.sd)[2]] <- paste0(cond, "_", time)
    }
  }
  # add the last condition with last time point
  count.matrix.av$new <- average_repl(count.matrix, "konost", 300)
  colnames(count.matrix.av)[dim(count.matrix.av)[2]] <- "konost_300"
  count.matrix.sd$new <- sd_repl(count.matrix, "konost", 300)
  colnames(count.matrix.sd)[dim(count.matrix.sd)[2]] <- "konost_300"
  return(list(count.matrix.av, count.matrix.sd))
}

create_plot_table <- function(count.matrix.av, count.matrix.sd) {
  # create table that contain time course data suitable for ggplot
  plot.data.av <- melt_data(count.matrix.av)
  plot.data.sd <- melt_data(count.matrix.sd)

  # merge average values with standard deviations
  plot.data.all <- merge(plot.data.av, plot.data.sd, by = c("id", "cond", "time"))
  colnames(plot.data.all) <- c("id", "cond", "time", "value", "sd")
  plot.data.all$ymax <- plot.data.all$value + plot.data.all$sd
  plot.data.all$ymin <- plot.data.all$value - plot.data.all$sd
  return(plot.data.all)
}

ave_repl <- function(count.matrix, name) {
  # change colnames of good replicates 3b_0, 3b_40, 4b_40 to 3_0, 3_40, 4_40
  colnames(count.matrix)[c(59, 65, 66)] <- c("wt_3_0", "wt_3_40", "wt_4_40")
  # define parameters
  conditions <- c("wt", "ko", "pten", "ki")

  res <- average_count_matrix(count.matrix, conditions)

  saveRDS(res[[1]],
    paste0("../pip3-rna-seq-output/rds/count-matrix-av", name, ".rds"))
  saveRDS(res[[2]],
    paste0("../pip3-rna-seq-output/rds/count-matrix-sd", name, ".rds"))

  plot.data.all <- create_plot_table(res[[1]], res[[2]])
  saveRDS(plot.data.all,
    paste0("../pip3-rna-seq-output/rds/plot-time-courses-all-genes", name, ".rds"))
}

pca <- function() {
  d <- readRDS("../pip3-rna-seq-output/rds/count-matrix-scaled.rds")
  res <- prcomp(d)

  pdf(file = "../pip3-rna-seq-output/figures/pca-variances.pdf", w=7, h=6)
  print(screeplot(res))
  dev.off()

  data <- as.data.frame(res$rotation[,1:3])
  data$cond <- unlist(lapply(strsplit(rownames(data), "_"), function(x){return(x[1])}))
  data$time <- unlist(lapply(strsplit(rownames(data), "_"), function(x){return(x[3])}))
  data$time <- as.numeric(data$time)
  data$cond <- factor(data$cond, levels = c("ko", "konost", "ki", "pten", "wt"))
  p <- ggplot(data, aes(PC1,PC2, color = cond)) +
      geom_point(aes(size = time, shape = cond)) +
      scale_size(range = c(3, 6)) +
      theme_bw()
  pdf(file = "../pip3-rna-seq-output/figures/pca12.pdf", w=7, h=6)
  print(p)
  dev.off()

  p <- ggplot(data, aes(PC2,PC3, color = cond)) +
      geom_point(aes(size = time, shape = cond)) +
      scale_size(range = c(3, 6)) +
      theme_bw()
  pdf(file = "../pip3-rna-seq-output/figures/pca23.pdf", w=7, h=6)
  print(p)
  dev.off()

  p <- ggplot(data, aes(PC1,PC3, color = cond)) +
      geom_point(aes(size = time, shape = cond)) +
      scale_size(range = c(3, 6)) +
      theme_bw()
  pdf(file = "../pip3-rna-seq-output/figures/pca13.pdf", w=7, h=6)
  print(p)
  dev.off()

  # in 3D
  open3d(windowRect=c(100,100,700,700))
  plot3d(res$rotation,xlab="PC1",ylab="PC2",zlab="PC3")
  spheres3d(res$rotation, radius=0.01,col=rainbow(length(res$rotation[,1])))
  grid3d(side="z", at=list(z=0))
  text3d(res$rotation, text=rownames(res$rotation), adj=1.3)
  rgl.postscript(file="../pip3-rna-seq-output/figures/pca-3d.pdf", fmt="pdf")
}
