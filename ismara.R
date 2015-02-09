process_averaged_tf_data <- function(file.name, sig.table) {
  d <- read.table(file = file.name, header = T)
  d$cond = rownames(d)
  d <- as.data.table(d)
  d[, time:=unlist(lapply(strsplit(d[,cond], "_"), "[[", 2))]
  d[, cond:=unlist(lapply(strsplit(d[,cond], "_"), "[[", 1))]
  d <- melt(d)
  d$time <- as.numeric(d$time)
  d$value <- as.numeric(d$value)
  d$group <- paste(d$cond, d$variable, sep = "_")
  d <- merge(d, sig.table, all = T)
  d$zval <- as.numeric(d$zval)
  d$variable <- factor(d$variable, levels = sig.table$variable)
  return(d)
}

process_tf_data <- function(file.name, sig.table) {

  v <- read.table(file = file.name, header = T)
  x <- data.frame(index = 1:75, names = rownames(v))
  x <- x[c(1,7,13,2,8,14,5,11,17,6,12,18,3,9,15,4,10,16,19,25,31,20,26,32,23,29,35,
    24,30,36,21,27,33,22,28,34,37,38,39,40,46,52,41,47,53,44,50,56,45,51,57,42,48,54,43,
    49,55,58,68,70,71,59,64,69,62,75,74,63,67,72,60,65,61,66,73),]
  v <- v[x$index, ]
  # v.split <- split(v, ceiling(seq_along(v[,1])/3))

  d <- data.frame()
  for(i in 1:length(v[1,])) {
    v.split <- split(v[,i], ceiling(seq_along(v[,i])/3))
    # this function is in svd.R
    t <- make_eigen_table_to_plot(v.split, colnames(v)[i])
    d <- rbind(d, t)
  }

  d <- merge(d, sig.table, all = T)

  d$time <- as.numeric(d$time)
  d$mean <- as.numeric(d$mean)
  d$sd <- as.numeric(d$sd)
  d$zval <- as.numeric(d$zval)
  d$group <- paste(d$cond, d$ind, sep = "_")
  d$ind <- factor(d$ind, levels = sig.table$ind)

  return(d)
}

plot_ismara_activities <- function() {
  # import table with TF motifs' significance values
  # note that I modified this file a bit before importing
  # (to match TF motif names in the activity table)
  sig <- read.table(file = "../pip3-rna-seq-input/ismara/active_matrices.txt")
  colnames(sig) <- c("ind", "zval")

  # import table with TF motifs' activities and standard deviations
  data <- process_tf_data("../pip3-rna-seq-input/ismara/activity_table.txt", sig)

  # calculate error bar values
  data$ymax <- data$mean + data$sd
  data$ymin <- data$mean - data$sd

  limits <- aes(ymax = ymax, ymin = ymin)

  # plot all TF motifs activities in one pdf file
  p <- ggplot(data, aes(time, mean, color = cond, group = group)) +
    geom_line() +
    facet_wrap( ~ ind, ncol = 10) +
    geom_errorbar(limits, width = 0.25)
  plots = dlply(data , "ind", `%+%`, e1 = p)
  ml = do.call(marrangeGrob, c(plots, list(nrow = 2, ncol = 1)))
  ggsave("../pip3-rna-seq-output/figures/ismara-activity-table.pdf", ml)
}

plot_ismara_averaged_activities <- function() {
  # import table with TF motifs' significance values
  # note that I modified this file a bit before importing
  # (to match TF motif names in the activity table)
  sig <- read.table(file = "../pip3-rna-seq-input/ismara/active_matrices_averaged.txt")
  colnames(sig) <- c("variable", "zval")

  # import table with TF motifs' activities and standard deviations
  data <- process_averaged_tf_data("../pip3-rna-seq-input/ismara/activity_table_averaged.txt", sig)
  err <- process_averaged_tf_data("../pip3-rna-seq-input/ismara/delta_table_averaged.txt", sig)

  # calculate error bar values
  data$ymax <- data$value + err$value
  data$ymin <- data$value - err$value

  limits <- aes(ymax = ymax, ymin = ymin)

  # plot all TF motifs activities in one pdf file
  p <- ggplot(data, aes(time, value, color = cond, group = group)) +
    geom_line() +
    facet_wrap( ~ variable, ncol = 10) +
    geom_errorbar(limits, width = 0.25)
  plots = dlply(data , "variable", `%+%`, e1 = p)
  ml = do.call(marrangeGrob, c(plots, list(nrow = 2, ncol = 1)))
  ggsave("../pip3-rna-seq-output/figures/ismara-activity-table-averaged.pdf", ml)
}

get_ismara_target_info_from_file <- function(file) {
  # read the targets file line by line and store them in mylist
  fc <- file(paste0("../pip3-rna-seq-input/ismara_targets/", file))
  mylist <- strsplit(readLines(fc), "\t")
  close(fc)

  # from mylist remove locations with no targets
  mylist <- mylist[unlist(lapply(mylist, length)) > 3]

  # extract information from mylist
  loc <- unlist(lapply(mylist, "[", 1))
  score <- unlist(lapply(mylist, "[", 2))
  motif <- unlist(lapply(mylist, "[", 3))
  len <- unlist(lapply(mylist, function(x) length(4:length(x))))
  targets <- unlist(lapply(mylist, function(x) x[4:length(x)]))
  target <- unlist(lapply(strsplit(targets, "\\|"), "[", 2))

  # put all information in one table
  mytable <- unique(data.table(location = rep(loc, len),
    score = rep(score, len), motif = rep(motif, len),
    target = target))
  return(mytable)
}

get_ismara_target_info_from_file_new <- function(file) {
  # read the targets file line by line and store them in mylist
  fc <- file(paste0("../pip3-rna-seq-input/ismara_targets/", file))
  mylist <- strsplit(readLines(fc), "\t")
  close(fc)

  # from mylist remove locations with no targets
  mylist <- mylist[unlist(lapply(mylist, length)) > 3]

  # extract information from mylist
  loc <- unlist(lapply(mylist, "[", 1))
  score <- unlist(lapply(mylist, "[", 2))
  motif <- unlist(lapply(mylist, "[", 3))
  len <- unlist(lapply(mylist, function(x) length(4:length(x))))
  targets <- unlist(lapply(mylist, function(x) x[4:length(x)]))
  target <- unlist(lapply(strsplit(targets, "\\|"), "[", 2))

  # record the source reference
  ref <- unlist(lapply(strsplit(targets, "\\|"), "[", 1))

  # put all information in one table
  mytable <- unique(data.table(location = rep(loc, len),
    score = rep(score, len), motif = rep(motif, len),
    target = target, ref = ref))

  # remove all GenBank entries, leave only RefSeq entries and collaps entries
  # with the same location,score,motif,target sets
  mytable <- unique(mytable[grepl("NM_", ref), list(location, score, motif, target)])
  return(mytable)
}


plot_motifs_tfs <- function(motif) {
  # read the targets file line by line and store them in mylist
  fc <- file(paste0("../pip3-rna-seq-input/ismara/motifs_tfs.txt"))
  mylist <- strsplit(readLines(fc), "\t")
  close(fc)

  entry <- mylist[[grep(motif, mylist)]]
  tfs <- entry[2:length(entry)]
  plot_genes(tfs, F, paste0("tfs-of-", motif))
}


merge_ismara_target_files <- function() {
  files <- list.files("../pip3-rna-seq-input/ismara_targets/")

  data <- data.table()
  for(i in files) {
    dt <- get_ismara_target_info_from_file(i)
    data <- rbindlist(list(data, dt))
  }

  data[, score:=as.numeric(data[,score])]

  # data <- data[, list(score = max(score)), by = c("target", "motif")]

  sig <- read.table(file = "../pip3-rna-seq-input/ismara/active_matrices.txt")
  colnames(sig) <- c("motif", "zval")

  data <- as.data.table(merge(as.data.frame(data), sig))

  saveRDS(data, "../pip3-rna-seq-output/rds/ismara-targets.rds")
}

merge_ismara_target_files_new <- function() {
  files <- list.files("../pip3-rna-seq-input/ismara_targets/")

  data <- data.table()
  for(i in files) {
    dt <- get_ismara_target_info_from_file_new(i)
    data <- rbindlist(list(data, dt))
  }

  data[, score:=as.numeric(data[,score])]

  data <- data[, list(score = max(score)), by = c("target", "motif")]

  sig <- read.table(file = "../pip3-rna-seq-input/ismara/active_matrices.txt")
  colnames(sig) <- c("motif", "zval")

  data <- as.data.table(merge(as.data.frame(data), sig))

  saveRDS(data, "../pip3-rna-seq-output/rds/ismara-targets-new.rds")
}

get_ismara_targets <- function(motifs) {
  d <- readRDS("../pip3-rna-seq-output/rds/ismara-targets.rds")
  d <- d[motif %in% motifs]
  target_ens <- hgnc_symbol_to_ensembl_id(d[,target])
  colnames(target_ens) <- c("target_ens", "target")
  target_ens <- as.data.table(target_ens)
  setkey(d, "target")
  setkey(target_ens, "target")
  # note that length res is < length d, because some gene names do not correspond
  # to hgnc names in my annotation table
  res <- d[target_ens]
  return(res)
}

get_ismara_motifs <- function(gene.set) {
  gene.names <- ensembl_id_to_hgnc_symbol(gene.set)
  d <- readRDS("../pip3-rna-seq-output/rds/ismara-targets.rds")
  return(d[target %in% gene.names$hgnc_symbol])
}

get_ismara_motifs_new <- function(gene.set) {
  gene.names <- ensembl_id_to_hgnc_symbol(gene.set)
  d <- readRDS("../pip3-rna-seq-output/rds/ismara-targets-new.rds")
  return(d[target %in% gene.names$hgnc_symbol])
}

ismara_plot_score_dist <- function(data, name) {
  ann.data <- data[1:6]
  ann.data$score <- 30
  ann.data$y <- Inf
  ann.data$target <- "target"
  ann.data$diff <- unique(ann.data$diff)[1]
  ann.data <- unique(ann.data)
  ann.data$D <- format(ann.data$D, digits = 1, nsmall = 2)

  # http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  cbPalette <- c("red", "blue")

  p <- ggplot(data, aes(x = score)) +
    geom_bar(aes(fill = diff), position = "dodge", width=.3) +
    # geom_density(aes(color = diff)) + 
    facet_wrap( ~ motif, scale = "free_y") +
    # geom_text(aes(x = score, y = y, label = D), data = ann.data, vjust = 1.2) +
    scale_fill_manual(values=cbPalette) +
    labs(x = "Target score, Sm", y = "Target number") +
    guides(fill = guide_legend(title = "Gene types")) +
    theme_bw()
    # theme(axis.line = element_blank(),
    #     panel.grid.major = element_blank(),
    #     panel.grid.minor = element_blank(),
    #     panel.border = element_blank(),
    #     panel.background = element_blank()) 
    # scale_x_log10()
  pdf(file = paste0("../pip3-rna-seq-output/figures/motif-diff-activity-", name, ".pdf"),
      w = 16, h = 10)
  print(p)
  dev.off()
}

sum_motif_target_scores_by_location <- function(t) {
  # sum all scores of each motif contributing to a given target
  t <- t[,list(score = sum(score)), by = c("motif","zval","target")]
  # filtering
  # t <- t[score > score.cut.off & zval > zval.cut.off]
  # calculate the number of target of each type for each motif and put it in one final
  # data.table d
  # t1 <- t[,list(diff.num = length(target)), by = c("motif","zval","diff")]
  # t2 <- t1[,list(tot.num = sum(diff.num)), by = c("motif","zval")]
  # t12 <- merge(t1, t2, by = c("motif", "zval"))
  # d <- merge(t, t12, by = c("motif", "zval", "diff"))
  return(t)
}

kolmogorov_smirnov_test <- function(score, diff, diff1, diff2) {
  x <- score[diff == diff1]
  # print(x)
  y <- score[diff == diff2]
  # print(y)
  if(length(x) != 0 & length(y) != 0) {
    return(ks.test(x, y))
  }
}

wilcox_test <- function(score, diff, diff1, diff2) {
  x <- score[diff == diff1]
  # print(x)
  y <- score[diff == diff2]
  # print(y)
  if(length(x) != 0 & length(y) != 0) {
    return(wilcox.test(x, y))
  }
}

motif_diff_activity <- function(gene.set1, gene.set2, diff1, diff2, name) {
  # get TF binding motifs (with zval > 0) from ISMARA report for act and pas genes
  d1 <- get_ismara_motifs(gene.set1)
  d1 <- sum_motif_target_scores_by_location(d1)

  d2 <- get_ismara_motifs(gene.set2)
  d2 <- sum_motif_target_scores_by_location(d2)

  d1 <- d1[,list(motif, zval, score, target, diff = diff1)]
  d2 <- d2[,list(motif, zval, score, target, diff = diff2)]

  d <- rbind(d1, d2)

  # filter insignificant targets and motifs (everywhere use 2)
  d <- d[zval > 2]
  
  ks <- d[,list(D = kolmogorov_smirnov_test(score, diff, diff1, diff2)[[1]],
                p = kolmogorov_smirnov_test(score, diff, diff1, diff2)[[2]]),
          by = "motif"]

  ks <- ks[order(-ks[,D])]
  # ks <- ks[p < 0.2]

  wx <- d[,list(W = wilcox_test(score, diff, diff1, diff2)[[1]],
                p = wilcox_test(score, diff, diff1, diff2)[[3]]),
          by = "motif"]

  wx <- wx[order(wx[,W])]
  # wx <- wx[p < 0.1]

  # all wx motifs are included in ks
  setkey(ks, "motif")
  setkey(d, "motif")
  d <- d[ks]

  d <- as.data.frame(d)
  # factorize motifs with levels proportional to their zval
  d$motif <- factor(d$motif, levels=unique(d[with(d, order(-D)), ]$motif))

  # plot distribution of scores for each motif
  ismara_plot_score_dist(d, name)
  return(as.data.table(d))
}

motif_diff_activity_new <- function(gene.set1, gene.set2, diff1, diff2, name) {
  # get TF binding motifs (with zval > 0) from ISMARA report for act and pas genes
  d1 <- get_ismara_motifs_new(gene.set1)
  # d1 <- sum_motif_target_scores_by_location(d1)

  d2 <- get_ismara_motifs_new(gene.set2)
  # d2 <- sum_motif_target_scores_by_location(d2)

  d1 <- d1[,list(motif, zval, score, target, diff = diff1)]
  d2 <- d2[,list(motif, zval, score, target, diff = diff2)]

  d <- rbind(d1, d2)

  # filter insignificant targets and motifs (everywhere use 2)
  d <- d[zval > 2]
  
  ks <- d[,list(D = kolmogorov_smirnov_test(score, diff, diff1, diff2)[[1]],
                p = kolmogorov_smirnov_test(score, diff, diff1, diff2)[[2]]),
          by = "motif"]

  ks <- ks[order(-ks[,D])]
  # ks <- ks[p < 0.2]

  wx <- d[,list(W = wilcox_test(score, diff, diff1, diff2)[[1]],
                p = wilcox_test(score, diff, diff1, diff2)[[3]]),
          by = "motif"]

  wx <- wx[order(wx[,W])]
  # wx <- wx[p < 0.1]

  # all wx motifs are included in ks
  setkey(ks, "motif")
  setkey(d, "motif")

  # use merge to keep motifs which present only in diff1 or diff2
  d <- merge(ks,d,all=TRUE, by = "motif")
  # motifs which present only in diff1 or diff2 will have score D=1
  d[is.na(D),D:=1]
  # d <- d[ks]

  d <- as.data.frame(d)
  # factorize motifs with levels proportional to their zval
  d$motif <- factor(d$motif, levels=unique(d[with(d, order(-D)), ]$motif))

  # plot distribution of scores for each motif
  ismara_plot_score_dist(d, name)
  return(as.data.table(d))
}

combine_gene_sets_motif_targets <- function(list, s, name) {
  i <- 1
  mod.list <- lapply(list, function(x) {
    tmp <- get_ismara_motifs_new(x)
    # tmp <- sum_motif_target_scores_by_location(tmp)
    tmp <- tmp[zval > 2 & score > s, list(motif, target, score)]
    tmp$set <- i
    tmp$x <- "1"
    i <<- i + 1
    setkey(tmp, "motif")
  })
  mod.list <- do.call("rbind", mod.list)
  mot.num <- length(unique(mod.list[,motif]))
  cols <- rainbow(mot.num, s=.6, v=.9)[sample(1:mot.num,mot.num)]
  p <- ggplot(mod.list, aes(x, x, color = motif)) +
    geom_jitter(aes(size = score), position = position_jitter(width = .2, height = 0)) +
    scale_size(range = c(2, 5)) +
    scale_colour_manual(values=cols) +
    # coord_flip() +
    facet_grid(target ~ set) +
    theme(strip.text.y = element_text(size = 7, angle = 0),
          panel.background = element_rect(fill = 'white', colour = 'gray'),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  pdf(file = paste0("../pip3-rna-seq-output/figures/motif-set-heatmap-", name, ".pdf"), w = 10, h = 15)
  print(p)
  dev.off()
  return(mod.list)
}

combine_motif_targets <- function(list, s, name) {
  i <- 1
  mod.list <- lapply(list, function(x) {
    tmp <- get_ismara_motifs(x)
    tmp <- sum_motif_target_scores_by_location(tmp)
    tmp <- tmp[zval > 2 & score > s, list(motif, target, score)]
    tmp$set <- i
    tmp$x <- "1"
    i <<- i + 1
    setkey(tmp, "motif")
  })
  mod.list <- do.call("rbind", mod.list)
  mot.num <- length(unique(mod.list[,motif]))
  cols <- rainbow(mot.num, s=.6, v=.9)[sample(1:mot.num,mot.num)]
  p <- ggplot(mod.list, aes(x, x, color = motif)) +
    geom_jitter(aes(size = score), position = position_jitter(width = .2, height = 0)) +
    scale_size(range = c(2, 5)) +
    scale_colour_manual(values=cols) +
    # coord_flip() +
    facet_grid(target ~ .) +
    theme(strip.text.y = element_text(size = 7, angle = 0),
          panel.background = element_rect(fill = 'white', colour = 'gray'),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  pdf(file = paste0("../pip3-rna-seq-output/figures/motif-set-heatmap-", name, ".pdf"), w = 6, h = 15)
  print(p)
  dev.off()
  return(mod.list)
}

split_targets_by_motifs <- function(targets, zv, s) {
  d <- get_ismara_motifs(targets)
  d <- d[score > s & zval > zv]
  motifs <- unique(d[, motif])
  for(m in motifs) {
    plot_motif_targets(d, m)
  }
}

plot_motif_targets <- function(target.table, m) {
  target.table <- sum_motif_target_scores_by_location(target.table)
  target.table <- target.table[motif == m, list(target, score)]
  # query biomart
  mart <- readRDS("../pip3-rna-seq-output/rds/count-matrix-ann.RDS")
  mart <- unique(mart[ , 1:2])
  colnames(mart)[2] <- "target"

  res <- merge(target.table, mart, by = "target")
  setnames(res, c("target", "score", "ensembl_gene_id"),
                c("name", "score", "id"))

  plot.data <- readRDS("../pip3-rna-seq-output/rds/plot-time-courses-all-genes.rds")
  plot.data <- merge(plot.data, res)

  plot.data$name <- factor(plot.data$name, 
    levels=unique(plot.data[with(plot.data, order(score)), ]$name))

  ann.data <- plot.data[, c(3:4, 8:9)]
  # ann.data$cond <- "wt"
  ann.data$time <- 200
  ann.data$value <- Inf
  ann.data <- unique(ann.data)
  ann.data$score <- format(ann.data$score, digits = 2, nsmall = 2)
  # ann.data$score <- as.numeric(ann.data$score)


  limits <- aes(ymax = ymax, ymin = ymin)

  p <- ggplot(plot.data, aes(time, value, group = cond, colour = cond)) +
      geom_line() +
      facet_wrap(~name, scale = "free_y") +
      geom_errorbar(limits, width = 0.25) +
      geom_text(aes(x = time, y = value, label = score),
        data = ann.data, vjust = 1.2, colour = "black")
  pdf(file = paste0("../pip3-rna-seq-output/figures/ismara-", m, "-targets.pdf"),
    width = (sqrt(length(res$name)) + 2)*3,
    height = sqrt(length(res$name))*3)
  print(p)
  dev.off()
}


plot_motif_diff_targets <- function(target.table, m) {
  target.table <- target.table[motif == m, list(diff, target, score)]
  # query biomart
  mart <- readRDS("../pip3-rna-seq-output/rds/count-matrix-ann.RDS")
  mart <- unique(mart[ , 1:2])
  colnames(mart)[2] <- "target"

  res <- merge(target.table, mart, by = "target")
  setnames(res, c("target", "diff", "score", "ensembl_gene_id"),
                c("name", "diff", "score", "id"))

  plot.data <- readRDS("../pip3-rna-seq-output/rds/plot-time-courses-all-genes.rds")
  plot.data <- merge(plot.data, res)

  plot.data$name <- factor(plot.data$name, 
    levels=unique(plot.data[with(plot.data, order(score)), ]$name))

  ann.data <- plot.data[, c(3:4, 8:10)]
  # ann.data$cond <- "wt"
  ann.data$time <- 200
  ann.data$value <- Inf
  ann.data <- unique(ann.data)
  ann.data$score <- format(ann.data$score, digits = 2, nsmall = 2)
  # ann.data$score <- as.numeric(ann.data$score)


  limits <- aes(ymax = ymax, ymin = ymin)

  p <- ggplot(plot.data,
    aes(time, value, group = cond, colour = cond)) +
      geom_line() + facet_wrap(name ~ diff, scale = "free_y") +
      geom_errorbar(limits, width = 0.25) +
      geom_text(aes(x = time, y = value, label = score),
        data = ann.data, vjust = 1.2, colour = "black")
  pdf(file = paste0("../pip3-rna-seq-output/figures/ismara-", m, "-targets.pdf"),
    width = (sqrt(length(res$name)) + 2)*3,
    height = sqrt(length(res$name))*3)
  print(p)
  dev.off()
}

plot_motif_target_hist <- function(targets, z, s) {
  d <- get_ismara_motifs(targets)
  d <- sum_motif_target_scores_by_location(d)
  d <- d[score > s & zval > z]
  t <- d[, list(target.num = length(target)), by = c("motif", "zval")]
  t <- t[order(zval)]
  # t <- t[zval > 2]
  # t <- t[order(zval)]
  t$motif <- factor(t$motif, levels = t$motif)
  p <- ggplot(t, aes(motif, target.num, fill = zval)) +
    geom_bar(stat = "identity") +
    coord_flip()
  return(list(data = d, plot = p))
}
