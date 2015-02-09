# update_GO_ann <- function() {
#   source("http://bioconductor.org/biocLite.R")
#   biocLite()
#   ann.BP <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "ensembl")
#   ann.MF <- annFUN.org("MF", mapping = "org.Hs.eg.db", ID = "ensembl")
#   ann.CC <- annFUN.org("CC", mapping = "org.Hs.eg.db", ID = "ensembl")

#   ann <- list(ann.BP, ann.MF, ann.CC)
#   t <- data.table()

#   for(j in ann) {
#     t <- rbind(t, rbindlist(lapply(seq_along(j), function(i) data.frame(gene.id = j[[i]], go.id = names(j)[[i]]))))
#   }
#   saveRDS(t, file = "../pip3-rna-seq-output/rds/GO-ann.rds")
# }

# g_profiler <- function(list.of.genes) {
#   file.remove("gprofiler.txt")
#   i <- 1
#   for (genes in list.of.genes) {
#     write(paste0("> eff", i, "\n", paste(genes, collapse = " ")),
#       file = "gprofiler.txt", append = T)
#     i <- i + 1;
#   }
# }

david <- function(list.of.genes) {
  write.table(as.data.frame(unlist(list.of.genes)), file = "david.txt",
    quote = F, row.names = F, col.names = F)
}

# create_go_data_frame <- function(name) {
#   files <- list.files("../pip3-rna-seq-output/GO")
#   files <- files[grepl(paste0(name, "-.-."), files)]

#   dat <- NULL

#   for(file in files) {
#     d <- tryCatch(read.table(paste0("../pip3-rna-seq-output/GO/", file, "/topgo-fisher-BP.txt")),
#      error=function(e) NULL)
#     # d <- read.table(paste0("../pip3-rna-seq-output/GO/", file, "/topgo-fisher-BP.txt"))
#     if(!is.null(d)) {
#       d[,3] <- file
#       dat <- rbind(dat, d)
#     }
#   }

#   colnames(dat) <- c("GO.Term.Accession", "pval", "set")
#   d <- read.csv("../pip3-rna-seq-input/annotations/gene-go-terms-GRCh37.p13.txt")
#   res <- merge(dat, d)
#   res[,2] <- as.numeric(res[,2])
#   saveRDS(res, paste0("../pip3-rna-seq-output/rds/go-terms-", name, ".rds"))
# }

# go_set_compare <- function(name, p.val) {
#   # create heatmaps of GO comparisons gene set clusters (pvalues inside)
#   res <- readRDS(paste0("../pip3-rna-seq-output/rds/go-terms-", name, ".rds"))
#   res <- res[res$pval < p.val,]

#   m <- unique(res[ , c(1,2,3)])
#   m <- dcast(m, GO.Term.Accession ~ set, value.var = "pval")

#   rownames(m) <- m$GO.Term.Accession
#   m <- m[,2:dim(m)[2]]
#   m <- as.matrix(m)
#   m[is.na(m)] <- 1

#   for(i in 1:dim(m)[2]) {
    
#     n <- colnames(m)[i]
#     d <- m[m[,i] != 1, ]
#     if(is.matrix(d)) {
#       d <- d[,colMeans(d) != 1]
#     }

#     if(is.matrix(d)) {
#       rownames(d) <- unique(res[res[,1] %in% rownames(d) , c(1,5)])[,2]
#       pdf(file = paste0("../pip3-rna-seq-output/figures/go-comp-", n, ".pdf"), w = 10, h = 12)
#       heatmap.2(log10(d), trace = "none", margins = c(6,22))
#       dev.off()
#     }
#   }
# }

# go_set <- function(name, set.ind, p.val) {
#   # create heatmaps of GO comparisons gene set clusters (pvalues inside)

#   res <- readRDS(paste0("../pip3-rna-seq-output/rds/go-terms-", name, ".rds"))
#   res <- res[res$pval < p.val,]

#   res <- res[grepl(paste0("-", set.ind, "-."), res$set), ]

#   m <- unique(res[ , c(2,3,5)])
#   m <- dcast(m, GO.Term.Name ~ set, value.var = "pval")

#   rownames(m) <- m$GO.Term.Name

#   if(dim(m)[2] > 2) {
#     m <- m[,2:dim(m)[2]]
#     m <- as.matrix(m)
#     m[is.na(m)] <- 1

#     pdf(file = paste0("../pip3-rna-seq-output/figures/go-", name, "-", set.ind, ".pdf"), w = 10, h = 25)
#     heatmap.2(log10(m), trace = "none", margins = c(6,22))
#     dev.off()
#   } else {
#     m <- m[order(m[,2]), ]
#     m[,1] <- factor(m[,1], levels = m[,1])
#     colnames(m)[2] <- "pval"
#     p <- ggplot(m, aes(GO.Term.Name, -log10(pval))) +
#       geom_bar(stat = "identity", width = 0.6) +
#       coord_flip() +
#       labs(x = "GO term name", y = "-log10(p-value)") +
#       theme(text = element_text(size = 16))
#     pdf(file = paste0("../pip3-rna-seq-output/figures/go-", name, "-", set.ind, ".pdf"), w = 10, h = 8)
#     print(p)
#     dev.off()
#   }
# }

# plot_go_terms <- function(name, stat, num) {
#   d <- read.table(paste0("../pip3-rna-seq-output/GO/", name, "/ann-BP-", stat, ".txt"), sep = "\t")
#   colnames(d) <- c("go.id", "go.name", "pval", "gene.id")
#   d <- as.data.table(d)

#   g.num <- d[,list(gene.num = length(unique(gene.id))), by = c("go.id", "pval",
#     "go.name")]
#   t <- tail(g.num[order(gene.num)], num)
#   t$go.name <- factor(t$go.name, levels = t$go.name)

#   p <- ggplot(t, aes(go.name, gene.num, fill = -log10(pval))) +
#     geom_bar(stat = "identity") +
#     coord_flip() +
#     # scale_fill_gradient2(midpoint = 10, low="blue", high="red") +
#     scale_fill_gradient2() +
#     theme_bw()

#   pdf(file = paste0("../pip3-rna-seq-output/figures/", name, "-BP-", stat, "-go-terms-by-num.pdf"),
#     w = 10, h = 0.18 * num)
#   print(p)
#   dev.off()

#   t <- head(g.num[order(pval)], num)
#   t$go.name <- factor(t$go.name, levels = rev(t$go.name))

#   p <- ggplot(t, aes(go.name, gene.num, fill = -log10(pval))) +
#     geom_bar(stat = "identity") +
#     coord_flip() +
#     scale_fill_gradient2() +
#     theme_bw()

#   pdf(file = paste0("../pip3-rna-seq-output/figures/", name, "-BP-", stat, "-go-terms-by-pval.pdf"),
#     w = 10, h = 0.18 * num)
#   print(p)
#   dev.off()

#   return(d)
# }

get_genes_from_go <- function(name, num) {
  d <- read.csv(paste0("../pip3-rna-seq-output/GO/", name, "/REVIGO.csv"))
  d <- d[,c(1, 2, 7)]
  d <- d[order(d$log10.p.value),]
  d <- as.data.table(d)
  go.terms <- d[1:num,term_ID]
  t <- readRDS(paste0("../pip3-rna-seq-output/GO/", name, "/genes-ann-BP.rds"))
  genes <- unique(unlist(t[names(t) %in% go.terms]))
  return(genes)
}

annotate_go_terms <- function(name) {
  d <- read.csv(paste0("../pip3-rna-seq-output/GO/", name, "/REVIGO.csv"))
  # d <- d[d$eliminated == 0, ]
  d <- d[,c(1, 2, 7)]
  d <- d[order(d$log10.p.value),]
  d <- as.data.table(d)
  setkey(d, "term_ID")

  t <- readRDS(paste0("../pip3-rna-seq-output/GO/", name, "/genes-ann-BP.rds"))
  t <- t[names(t) %in% d$term_ID]
  g.num <- lapply(t, length)
  g.num <- as.data.frame(unlist(g.num))
  g.num$term_ID <- rownames(g.num)
  colnames(g.num)[1] <- "num"
  g.num <- as.data.table(g.num)
  setkey(g.num, "term_ID")
  res <- d[g.num]
  return(res)
}

plot_go_terms <- function(d, name, num) {
  d <- d[order(log10.p.value)]
  d$description <- factor(d$description, levels = rev(d$description))

  d <- head(d, num)

  p <- ggplot(d, aes(description, num, fill = log10.p.value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    # scale_fill_gradient2(midpoint = 10, low="blue", high="red") +
    scale_fill_gradient(high = "grey90", low = "grey40") +
    theme_bw() +
    theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
  pdf(file = paste0("../pip3-rna-seq-output/figures/go-", name, ".pdf"),
    w = 10, h = 4)
  print(p)
  dev.off()
}

# plot_go_terms2 <- function(d1, d2, name, num) {
#   d1 <- d1[order(log10.p.value)]
#   d1$description <- factor(d1$description, levels = rev(d1$description))
#   d1 <- head(d1, num)
#   d1$name <- "8-1"

#   d2 <- d2[order(log10.p.value)]
#   d2$description <- factor(d2$description, levels = rev(d2$description))
#   d2 <- head(d2, num)
#   d2$name <- "8-2"

#   d <- rbind(d1, d2)

#   p <- ggplot(d, aes(description, num, fill = log10.p.value)) +
#     geom_bar(stat = "identity") +
#     coord_flip() +
#     # scale_fill_gradient2(midpoint = 10, low="blue", high="red") +
#     scale_fill_gradient2() +
#     theme_bw() +
#     facet_grid(. ~ name, scale = "free_y")
#   pdf(file = paste0("../pip3-rna-seq-output/figures/go-", name, ".pdf"),
#     w = 10, h = 4)
#   print(p)
#   dev.off()
# }

quickgo1 <- function(go.set, name, num) {
  t <- head(go.set[order(log10.p.value)], num)

  t$color <- colorRampPalette(c("blue", "lightblue"))(num)

  system(paste0("rm -r ../pip3-rna-seq-output/quickgo/", name))
  system(paste0("mkdir ../pip3-rna-seq-output/quickgo/", name))

  pdf(file = paste0("../pip3-rna-seq-output/quickgo/", name, "/bl-grad.pdf"))
  plot(rep(1,num),col=colorRampPalette(c("blue", "lightblue"))(num), pch=19, cex=3)
  dev.off()

  sink(paste0("../pip3-rna-seq-output/quickgo/", name, "/goTerms.txt"))
  for(i in 1:dim(t)[1]) {
    cat(paste0(t[i]$color, "\n"))
    cat(paste0(t[i]$term_ID, "\n"))
  }
  sink()

  setwd("quickGO5.1-BI")
  system(paste0("sh createGoChart.sh ../../pip3-rna-seq-output/quickgo/", name, "/goTerms.txt"))
  system(paste0("cp goGraphProcess.png ../../pip3-rna-seq-output/quickgo/", name))
  system(paste0("cp goGraphProcess.svg ../../pip3-rna-seq-output/quickgo/", name))
  setwd("..")
  # return(t)
}

quickgo2 <- function(go.set1, go.set2, name, num) {
  t1 <- head(go.set1[order(log10.p.value)], num)
  t2 <- head(go.set2[order(log10.p.value)], num)
  # bluefunc <- colorRampPalette(c("blue", "white"))
  # redfunc <- colorRampPalette(c("red", "white"))

  # bl = bluefunc(num)[findInterval(-t1$log10.p.value, seq(1:num))]
  # rd = redfunc(num)[findInterval(-t2$log10.p.value, seq(1:num))]

  t1$color <- colorRampPalette(c("blue", "lightblue"))(num)
  t2$color <- colorRampPalette(c("darkgreen", "green"))(num)

  system(paste0("rm -r ../pip3-rna-seq-output/quickgo/", name))
  system(paste0("mkdir ../pip3-rna-seq-output/quickgo/", name))

  pdf(file = paste0("../pip3-rna-seq-output/quickgo/", name, "/bl-grad.pdf"))
  plot(rep(1,num),col=colorRampPalette(c("blue", "lightblue"))(num), pch=19, cex=3)
  dev.off()
  pdf(file = paste0("../pip3-rna-seq-output/quickgo/", name, "/gr-grad.pdf"))
  plot(rep(1,num),col=colorRampPalette(c("darkgreen", "green"))(num), pch=19, cex=3)
  dev.off()

  t1.out <- t1[,list(term_ID, log10.p.value, color)]
  t2.out <- t2[,list(term_ID, log10.p.value, color)]

  setkey(t1.out, "term_ID")
  setkey(t2.out, "term_ID")

  t.out <- merge(t1.out, t2.out, all = T)
  t.out[is.na(color.x), color:=color.y]
  t.out[is.na(color.y), color:=color.x]
  t.out[is.na(color), color:="#FE2E2E"]

  sink(paste0("../pip3-rna-seq-output/quickgo/", name, "/goTerms.txt"))
  for(i in 1:dim(t.out)[1]) {
    cat(paste0(t.out[i]$color, "\n"))
    cat(paste0(t.out[i]$term_ID, "\n"))
  }
  sink()

  setwd("quickGO5.1-BI")
  system(paste0("sh createGoChart.sh ../../pip3-rna-seq-output/quickgo/", name, "/goTerms.txt"))
  system(paste0("cp goGraphProcess.png ../../pip3-rna-seq-output/quickgo/", name))
  system(paste0("cp goGraphProcess.svg ../../pip3-rna-seq-output/quickgo/", name))
  setwd("..")
  # return(t.out)
}

GO <- function(selected.genes, all.genes, name, p.val) {
  system(paste0("rm -r ../pip3-rna-seq-output/GO/", name))
  system(paste0("mkdir ../pip3-rna-seq-output/GO/", name))

  geneList <- factor(as.integer(all.genes %in% selected.genes))
  names(geneList) <- all.genes

  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
    annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "ensembl")

  res <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

  res.table <- GenTable(GOdata, classicFisher = res,
    topNodes = length(usedGO(GOdata)))
  res.table$classicFisher <- as.numeric(res.table$classicFisher)

  res.table <- res.table[res.table$classicFisher < p.val,]

  genes.ann <- genesInTerm(GOdata, whichGO=res.table$GO.ID)
  genes.ann <- lapply(genes.ann, function(x) {x[x %in% selected.genes]})
  saveRDS(genes.ann, paste0("../pip3-rna-seq-output/GO/", name, "/genes-ann-BP.rds"))

  # write.table(res.table, file = paste0("../pip3-rna-seq-output/GO/", name, "/BP-Fisher.txt"),
  #   quote = F, row.names = F, col.names = F, sep = "\t")

  write.table(res.table[,c(1,6)], file = paste0("../pip3-rna-seq-output/GO/", name, "/BP-Fisher-revigo.txt"),
    quote = F, row.names = F, col.names = F, sep = "\t")

  # printGraph(GOdata, res, firstSigNodes = 20,
  #   fn.prefix = paste0("../pip3-rna-seq-output/GO/", name, "/Fisher"),
  #   useInfo = "def", pdfSW = TRUE)

  system(paste0("python2.7 revigo.py BP-Fisher-revigo.txt ", name))
}
