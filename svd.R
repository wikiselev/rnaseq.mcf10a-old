############################
# functions for SVD analysis
############################

# ----- Define a function for plotting a matrix ----- #
# taken from here http://www.phaget4.org/R/image_matrix.html
myImagePlot <- function(x, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
       yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
       xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
       title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
     xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
     yLabels <- c(1:nrow(x))
  }

  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))

  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0,1,length=256),  # Red
                   seq(0,1,length=256),  # Green
                   seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))

  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]

  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
  ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
  cex.axis=0.7)

  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n")

  layout(1)
}
# ----- END plot function ----- #

plot_sample_data <- function(sample.size, m, d, u, v) {
  genes <- sample(rownames(m), sample.size)
  m <- m[rownames(m) %in% genes, ]
  u <- u[rownames(m) %in% genes, ]

  pdf(file = "../pip3-rna-seq-output/figures/sm.pdf")
  myImagePlot(m)
  dev.off()

  pdf(file = "../pip3-rna-seq-output/figures/u.pdf")
  myImagePlot(u)
  dev.off()

  pdf(file = "../pip3-rna-seq-output/figures/v.pdf")
  myImagePlot(v)
  dev.off()

  pdf(file = "../pip3-rna-seq-output/figures/d.pdf")
  myImagePlot(diag(d))
  dev.off()
}

make_eigen_table_to_plot <- function(v.split, v.ind) {
  d <- data.frame(ind = v.ind, cond = "PI3K KI", time = 0, mean = mean(v.split$`1`),
    sd = sd(v.split$`1`), stringsAsFactors=FALSE)
  d <- rbind(d, c(v.ind, "PI3K KI", 15, mean(v.split$`2`), sd(v.split$`2`)))
  d <- rbind(d, c(v.ind, "PI3K KI", 40, mean(v.split$`3`), sd(v.split$`3`)))
  d <- rbind(d, c(v.ind, "PI3K KI", 90, mean(v.split$`4`), sd(v.split$`4`)))
  d <- rbind(d, c(v.ind, "PI3K KI", 180, mean(v.split$`5`), sd(v.split$`5`)))
  d <- rbind(d, c(v.ind, "PI3K KI", 300, mean(v.split$`6`), sd(v.split$`6`)))
  d <- rbind(d, c(v.ind, "A66", 0, mean(v.split$`7`), sd(v.split$`7`)))
  d <- rbind(d, c(v.ind, "A66", 15, mean(v.split$`8`), sd(v.split$`8`)))
  d <- rbind(d, c(v.ind, "A66", 40, mean(v.split$`9`), sd(v.split$`9`)))
  d <- rbind(d, c(v.ind, "A66", 90, mean(v.split$`10`), sd(v.split$`10`)))
  d <- rbind(d, c(v.ind, "A66", 180, mean(v.split$`11`), sd(v.split$`11`)))
  d <- rbind(d, c(v.ind, "A66", 300, mean(v.split$`12`), sd(v.split$`12`)))
  d <- rbind(d, c(v.ind, "A66 (no EGF)", 300, mean(v.split$`13`), sd(v.split$`13`)))
  d <- rbind(d, c(v.ind, "PTEN KO", 0, mean(v.split$`14`), sd(v.split$`14`)))
  d <- rbind(d, c(v.ind, "PTEN KO", 15, mean(v.split$`15`), sd(v.split$`15`)))
  d <- rbind(d, c(v.ind, "PTEN KO", 40, mean(v.split$`16`), sd(v.split$`16`)))
  d <- rbind(d, c(v.ind, "PTEN KO", 90, mean(v.split$`17`), sd(v.split$`17`)))
  d <- rbind(d, c(v.ind, "PTEN KO", 180, mean(v.split$`18`), sd(v.split$`18`)))
  d <- rbind(d, c(v.ind, "PTEN KO", 300, mean(v.split$`19`), sd(v.split$`19`)))
  d <- rbind(d, c(v.ind, "WT", 0, mean(v.split$`20`), sd(v.split$`20`)))
  d <- rbind(d, c(v.ind, "WT", 15, mean(v.split$`21`), sd(v.split$`21`)))
  d <- rbind(d, c(v.ind, "WT", 40, mean(v.split$`22`), sd(v.split$`22`)))
  d <- rbind(d, c(v.ind, "WT", 90, mean(v.split$`23`), sd(v.split$`23`)))
  d <- rbind(d, c(v.ind, "WT", 180, mean(v.split$`24`), sd(v.split$`24`)))
  d <- rbind(d, c(v.ind, "WT", 300, mean(v.split$`25`), sd(v.split$`25`)))
  return(d)
}

plot_eigen_arrays_6_points <- function(v, name) {
  d <- data.frame()
  for(i in 1:5) {
    v.split <- split(v[,i], ceiling(seq_along(v[,i])/3))
    t <- make_eigen_table_to_plot(v.split, paste0("V", i))
    d <- rbind(d, t)
  }

  d$mean <- as.numeric(d$mean)
  d$sd <- as.numeric(d$sd)
  d$time <- as.numeric(d$time)

  limits <- aes(ymax = mean + sd, ymin = mean - sd)

  p <- ggplot(d, aes(time, mean, colour = cond, group = cond)) + geom_line() +
    geom_errorbar(limits) +
    facet_wrap( ~ ind, scale = "free") +
    labs(x = "Time, min", y = "") +
    theme(legend.title=element_blank())

  pdf(file = paste0("../pip3-rna-seq-output/figures/", name, ".pdf"), w = 12, h = 6)
  print(p)
  dev.off()
}

plot_eigen_arrays_4_points <- function(v, name) {
  v.split <- split(v, ceiling(seq_along(v)/3))
  d <- data.frame(cond = "ki", time = 0, mean = mean(v.split$`1`),
    sd = sd(v.split$`1`), stringsAsFactors=FALSE)
  d <- rbind(d, c("ki", 15, mean(v.split$`2`), sd(v.split$`2`)))
  d <- rbind(d, c("ki", 40, mean(v.split$`3`), sd(v.split$`3`)))
  d <- rbind(d, c("ki", 90, mean(v.split$`4`), sd(v.split$`4`)))
  d <- rbind(d, c("a66", 0, mean(v.split$`5`), sd(v.split$`5`)))
  d <- rbind(d, c("a66", 15, mean(v.split$`6`), sd(v.split$`6`)))
  d <- rbind(d, c("a66", 40, mean(v.split$`7`), sd(v.split$`7`)))
  d <- rbind(d, c("a66", 90, mean(v.split$`8`), sd(v.split$`8`)))
  d <- rbind(d, c("pten", 0, mean(v.split$`9`), sd(v.split$`9`)))
  d <- rbind(d, c("pten", 15, mean(v.split$`10`), sd(v.split$`10`)))
  d <- rbind(d, c("pten", 40, mean(v.split$`11`), sd(v.split$`11`)))
  d <- rbind(d, c("pten", 90, mean(v.split$`12`), sd(v.split$`12`)))
  d <- rbind(d, c("wt", 0, mean(v.split$`13`), sd(v.split$`13`)))
  d <- rbind(d, c("wt", 15, mean(v.split$`14`), sd(v.split$`14`)))
  d <- rbind(d, c("wt", 40, mean(v.split$`15`), sd(v.split$`15`)))
  d <- rbind(d, c("wt", 90, mean(v.split$`16`), sd(v.split$`16`)))

  d$mean <- as.numeric(d$mean)
  d$sd <- as.numeric(d$sd)
  d$time <- as.numeric(d$time)

  limits <- aes(ymax = mean + sd, ymin = mean - sd)

  p <- ggplot(d, aes(time, mean, colour = cond, group = cond)) + geom_line() +
    geom_errorbar(limits)
  pdf(file = paste0("../pip3-rna-seq-output/figures/", name, ".pdf"))
  print(p)
  dev.off()
}

plot_projection_all <- function(m, v) {
  # scalar multiplication is a projection of a vector to another vector
  proj <- as.data.frame(m%*%v)
  proj1 <- melt(proj)
  p <- ggplot(proj1, aes(value, color = variable)) + geom_density()
  pdf(file = paste0("../pip3-rna-seq-output/figures/projection-all.pdf"))
  print(p)
  dev.off()
  return(proj)
}

plot_variance <- function(d) {
  d <- data.frame(x = paste("V", c(1:length(d)), sep = ""),
    y = d**2/sum(d**2)*100)
  d <- d[1:5,]
  p <- ggplot(d, aes(as.factor(x), y)) +
    geom_bar(stat = "identity", width = 0.4) +
    labs(x = "Right singular vectors", y = "% of data variance")
  pdf(file = "../pip3-rna-seq-output/figures/svd-variance.pdf", width = 6, height = 4)
  print(p)
  dev.off()
}

split_genes <- function(m) {
  # retrieving annotations from ensembl
  annotations <- data.frame()
  # biomaRt attributes
  attr <- c("ensembl_gene_id", 
        "gene_biotype")
  annotations <- getBM(attributes = attr, mart = mart)
  prot.coding <-
    m[rownames(m) %in%
      annotations[annotations$gene_biotype == "protein_coding", 1], ]
  other <-
    m[!(rownames(m) %in%
      annotations[annotations$gene_biotype == "protein_coding", 1]), ]
  return(list(prot.coding = prot.coding, other = other))
}

s_v_d <- function(m, num.of.points) {
  # make a singular value decomposition
  m.svd <- svd(m)
  # names(m.svd)
  # [1] "d" "u" "v"
  # d - singularity diagonal
  # v - eigentarrays
  # u - eigengenes

  d <- m.svd$d
  # length(d)
  # [1] 75
  u <- m.svd$u
  # dim(u)
  # [1] 38662 75
  v <- m.svd$v
  # dim(v)
  # [1] 75  75

  # plot_sample_data(100, m, d, u, v)
  plot_variance(d)

  if(num.of.points == 6) {
    plot_eigen_arrays_6_points(v, "svd-right-singular-vectors")
    # plot_eigen_arrays_6_points(v[,2], "eigen_array2")
    # plot_eigen_arrays_6_points(v[,3], "eigen_array3")
    # plot_eigen_arrays_6_points(v[,4], "eigen_array4")
    # plot_eigen_arrays_6_points(v[,5], "eigen_array5")
  } else {
    plot_eigen_arrays_4_points(v[,1], "eigen-array1")
    plot_eigen_arrays_4_points(v[,2], "eigen-array2")
    plot_eigen_arrays_4_points(v[,3], "eigen-array3")
    plot_eigen_arrays_4_points(v[,4], "eigen-array4")
    plot_eigen_arrays_4_points(v[,5], "eigen-array5")
  }
  return(v[,1:5])
}

lemontree_make_expr_matrix <- function(genes) {
  # import count matrix
  cm <- readRDS("../pip3-rna-seq-output/rds/count-matrix.rds")
  # select genes of interest
  cm <- cm[rownames(cm) %in% genes,]
  # scale count matrix
  cm <- t(scale(t(cm), center = T, scale = T))
  write.table(cm, file = "../lemontree-test/expr_matrix.txt", sep = "\t", quote = F)
}

s_v_d_genes <- function(genes, num.of.points) {
  # import count matrix
  cm <- readRDS("../pip3-rna-seq-output/rds/count-matrix.rds")

  if(num.of.points == 4) {
    cm <- cm[, !grepl("_300", colnames(cm))]
    cm <- cm[, !grepl("_180", colnames(cm))]
  }
  # select genes of interest
  cm <- cm[rownames(cm) %in% genes,]
  # scale count matrix
  cm <- t(scale(t(cm), center = T, scale = T))
  # perform an svd analysis on count.matrix
  eigenarr <- s_v_d(cm, num.of.points)
  # project data to new eigenarrays
  res <- cm %*% eigenarr
  return(res)
}

s_v_d_ismara <- function() {
  # svd analysis of raw TF-activity data
  d <- t(read.table(file = "../pip3-rna-seq-input/ismara/activity_table.txt", header = T))
  # scale count matrix
  d <- t(scale(t(d), center = T, scale = T))
  # perform an svd analysis on activity matrix
  eigenarr <- s_v_d(d, 6)
  # project data to new eigenarrays
  res <- d %*% eigenarr

  sig <- read.table(file = "../pip3-rna-seq-input/ismara/active_matrices.txt")
  colnames(sig) <- c("variable", "zval")
  # zval of 1.5 was taken from Luisier et al., 2014
  sig.z <- sig[sig$zval > 2.0,]
  res <- res[rownames(res) %in% sig.z[,1],]

  pdf(file = "../pip3-rna-seq-output/figures/tf-heatmap.pdf", w = 10, h = 8)
  heatmap.3(abs(res), Rowv = T, cluster.by.row = T, Colv = F,
    cluster.by.col = F, dendrogram = "row", breaks = 6,
    main = "Absolute value of TF binding motif activity projections",
    font.main = 1)
  dev.off()
  return(res)
}

arrange_ismara_activity_matrix <- function(count.matrix){
        x <- c(37,38,39,
               58,68,70,59,64,71,62,69,75,63,67,74,60,65,72,61,66,73,
               19,25,31,20,26,32,23,29,35,24,30,36,21,27,33,22,28,34,
               1,7,13,2,8,14,5,11,17,6,12,18,3,9,15,4,10,16,
               40,46,52,41,47,53,44,50,56,45,51,57,42,48,54,43,49,55)
        return(count.matrix[,x])
}

heatmap_ismara_activities <- function() {
        d <- t(read.table(file = "../pip3-rna-seq-input/ismara/activity_table.txt", header = T))
        sig <- read.table(file = "../pip3-rna-seq-input/ismara/active_matrices.txt")
        colnames(sig) <- c("variable", "zval")
        sig.z <- sig[sig$zval > 2.0,]
        d <- d[rownames(d) %in% sig.z[,1],]
        t <- arrange_ismara_activity_matrix(d)
        heatmap.2(t, Colv=F, col=redgreen(19), trace = "none", dendrogram = "row")
}
