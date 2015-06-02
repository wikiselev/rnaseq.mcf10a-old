annotate_go_terms <- function(name) {
    d <- read.csv(paste0("../pip3-rna-seq-output/GO/", name, "/REVIGO.csv"))
    # d <- d[d$eliminated == 0, ]
    d <- d[, c(1, 2, 7)]
    d <- d[order(d$log10.p.value), ]
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
    
    p <- ggplot(d, aes(description, num, fill = log10.p.value)) + geom_bar(stat = "identity") + coord_flip() + 
        # scale_fill_gradient2(midpoint = 10, low='blue', high='red') +
    scale_fill_gradient(high = "grey90", low = "grey40") + theme_bw() + theme(axis.line = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        panel.background = element_blank())
    pdf(file = paste0("../pip3-rna-seq-output/figures/go-", name, ".pdf"), w = 10, h = 4)
    print(p)
    dev.off()
}

GO <- function(selected.genes, all.genes, name, p.val) {
    system(paste0("rm -r ../pip3-rna-seq-output/GO/", name))
    system(paste0("mkdir ../pip3-rna-seq-output/GO/", name))
    
    geneList <- factor(as.integer(all.genes %in% selected.genes))
    names(geneList) <- all.genes
    
    GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.org, mapping = "org.Hs.eg.db", 
        ID = "ensembl")
    
    res <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    
    res.table <- GenTable(GOdata, classicFisher = res, topNodes = length(usedGO(GOdata)))
    res.table$classicFisher <- as.numeric(res.table$classicFisher)
    
    res.table <- res.table[res.table$classicFisher < p.val, ]
    
    genes.ann <- genesInTerm(GOdata, whichGO = res.table$GO.ID)
    genes.ann <- lapply(genes.ann, function(x) {
        x[x %in% selected.genes]
    })
    saveRDS(genes.ann, paste0("../pip3-rna-seq-output/GO/", name, "/genes-ann-BP.rds"))
    
    # write.table(res.table, file = paste0('../pip3-rna-seq-output/GO/', name, '/BP-Fisher.txt'), quote = F,
    # row.names = F, col.names = F, sep = '\t')
    
    write.table(res.table[, c(1, 6)], file = paste0("../pip3-rna-seq-output/GO/", name, "/BP-Fisher-revigo.txt"), 
        quote = F, row.names = F, col.names = F, sep = "\t")
    
    # printGraph(GOdata, res, firstSigNodes = 20, fn.prefix = paste0('../pip3-rna-seq-output/GO/', name,
    # '/Fisher'), useInfo = 'def', pdfSW = TRUE)
    
    system(paste0("python2.7 revigo.py BP-Fisher-revigo.txt ", name))
} 
