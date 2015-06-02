# hier_clust <- function(genes) { cm <- readRDS('../pip3-rna-seq-output/rds/count-matrix-scaled.rds') #
# select genes of interest cm <- cm[rownames(cm) %in% genes,] h.c <- hclust(dist(cm), method = 'complete')
# }

# k_means_clust <- function(genes, clust.num, name) { set.seed(2) cm <-
# readRDS('../pip3-rna-seq-output/rds/count-matrix-scaled.rds') # select genes of interest cm <-
# cm[rownames(cm) %in% genes,] k.m <- kmeans(cm, clust.num, nstart = 20) plot_clust_profile(k.m$cluster,
# paste0(name, '-k-m-', clust.num)) saveRDS(k.m, paste0('../pip3-rna-seq-output/rds/', name, '-k-m-',
# clust.num, '.rds')) # return(k.m) }

clust_boot <- function(genes, min.clust.num, max.clust.num, name, ind) {
    cm <- readRDS("../pip3-rna-seq-output/rds/count-matrix-scaled.rds")
    # select genes of interest
    cm <- cm[rownames(cm) %in% genes, ]
    for (i in min.clust.num:max.clust.num) {
        c.boot <- clusterboot(cm, B = 100, bootmethod = "boot", clustermethod = pamkCBI, krange = i, seed = 155)
        saveRDS(c.boot, paste0("../pip3-rna-seq-output/rds/clusts/", name, "/", ind, "-", i, ".rds"))
    }
}

clust_boot_a66 <- function(genes, min.clust.num, max.clust.num, name, ind) {
    cm <- readRDS("../pip3-rna-seq-output/rds/count-matrix-scaled.rds")
    
    cm <- cm[, c(19:39, 58:75)]
    # select genes of interest
    cm <- cm[rownames(cm) %in% genes, ]
    for (i in min.clust.num:max.clust.num) {
        c.boot <- clusterboot(cm, B = 100, bootmethod = "boot", clustermethod = pamkCBI, krange = i, seed = 155)
        saveRDS(c.boot, paste0("../pip3-rna-seq-output/rds/clusts/", name, "/", ind, "-", i, ".rds"))
    }
}


plot_bootstrap_data <- function(name) {
    files <- list.files(paste0("../pip3-rna-seq-output/rds/clusts/", name))
    # files <- files[grepl(paste0('clust-', name, 'set'), files)]
    jaccard.stats <- data.frame(factor(), numeric(), numeric(), factor(), factor())
    for (file in files) {
        set <- strsplit(file, "-")[[1]][1]
        clust.num <- as.numeric(strsplit(strsplit(file, "-")[[1]][2], "\\.")[[1]][1])
        dat <- readRDS(paste0("../pip3-rna-seq-output/rds/clusts/", name, "/", file))
        for (i in 1:clust.num) {
            jaccard.stats <- rbind(jaccard.stats, cbind(set, dat$bootresult[i, ], i, clust.num))
        }
        
    }
    colnames(jaccard.stats) <- c("set", "jacc.sim", "clust.ind", "clust.num")
    jaccard.stats$jacc.sim <- as.numeric(jaccard.stats$jacc.sim)
    
    # limits <- aes(ymax = jacc.sim + jacc.sd, ymin = jacc.sim - jacc.sd) jaccard.stats$set <-
    # factor(jaccard.stats$set, levels = unique(jaccard.stats$set))
    
    p <- ggplot(jaccard.stats, aes(clust.ind, jacc.sim)) + geom_boxplot() + facet_grid(set ~ clust.num) + 
        geom_hline(yintercept = 0.75, color = "red") + theme_bw()
    pdf(file = paste0("../pip3-rna-seq-output/figures/cluster-", name, "-boots.pdf"), width = 12, height = 8)
    print(p)
    dev.off()
}

plot_all_clusts <- function(clusts, name, inds) {
    plot.data <- readRDS("../pip3-rna-seq-output/rds/plot-time-courses-all-genes-scaled.rds")
    
    plot.data <- as.data.table(plot.data)
    plot.data[, `:=`(clust, -1)]
    
    ind <- 1
    clust.size <- NULL
    
    for (j in clusts) {
        clust.size <- rbind(clust.size, cbind(as.data.frame(table(j)), inds[ind]))
        for (i in 1:length(unique(j))) {
            genes <- names(j[j == i])
            plot.data[id %in% genes, `:=`(clust, i)]
            plot.data[id %in% genes, `:=`(set, inds[ind])]
        }
        ind <- ind + 1
    }
    
    colnames(clust.size) <- c("clust", "gene.num", "set")
    clust.size$val <- 0
    clust.size$sd <- 0
    clust.size$cond <- "wt"
    clust.size$time <- 0
    
    plot.data <- plot.data[clust != -1]
    # plot.data[,gr:=paste(id, cond, sep = '_')]
    
    plot.data <- plot.data[, list(val = mean(value), sd = sd(value)), by = c("cond", "time", "clust", "set")]
    
    limits <- aes(ymax = val + sd, ymin = val - sd)
    
    p <- ggplot(plot.data, aes(time, val, group = cond, color = cond)) + geom_line() + geom_point() + geom_text(aes(x = 270, 
        y = 1, label = gene.num), data = clust.size, colour = "black") + facet_grid(set ~ clust) + theme_bw()
    
    pdf(file = paste0("../pip3-rna-seq-output/figures/clusts-", name, "-av-all.pdf"), width = 10, height = 10)
    print(p)
    dev.off()
}

plot_all_clusts_a66 <- function(clusts, name, inds) {
    plot.data <- readRDS("../pip3-rna-seq-output/rds/plot-time-courses-all-genes-scaled.rds")
    
    plot.data <- as.data.table(plot.data)
    
    plot.data <- plot.data[!cond %in% c("pten", "ki")]
    
    plot.data[, `:=`(clust, -1)]
    
    ind <- 1
    clust.size <- NULL
    
    for (j in clusts) {
        clust.size <- rbind(clust.size, cbind(as.data.frame(table(j)), inds[ind]))
        for (i in 1:length(unique(j))) {
            genes <- names(j[j == i])
            plot.data[id %in% genes, `:=`(clust, i)]
            plot.data[id %in% genes, `:=`(set, inds[ind])]
        }
        ind <- ind + 1
    }
    
    colnames(clust.size) <- c("clust", "gene.num", "set")
    clust.size$val <- 0
    clust.size$sd <- 0
    clust.size$cond <- "wt"
    clust.size$time <- 0
    
    plot.data <- plot.data[clust != -1]
    # plot.data[,gr:=paste(id, cond, sep = '_')]
    
    plot.data <- plot.data[, list(val = mean(value), sd = sd(value)), by = c("cond", "time", "clust", "set")]
    
    limits <- aes(ymax = val + sd, ymin = val - sd)
    
    p <- ggplot(plot.data, aes(time, val, group = cond, color = cond)) + geom_line() + geom_point() + geom_text(aes(x = 270, 
        y = 1, label = gene.num), data = clust.size, colour = "black") + facet_grid(set ~ clust, scale = "free_y")
    
    pdf(file = paste0("../pip3-rna-seq-output/figures/clusts-", name, "-av-all.pdf"), width = 10, height = 10)
    print(p)
    dev.off()
}

get_clust_genes <- function(name, file.name) {
    t <- readRDS(paste0("../pip3-rna-seq-output/rds/clusts/", name, "/", file.name, ".rds"))
    return(t)
} 
