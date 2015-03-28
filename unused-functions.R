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

plot_motifs_tfs <- function(motif) {
        # read the targets file line by line and store them in mylist
        fc <- file(paste0("../pip3-rna-seq-input/ismara/motifs_tfs.txt"))
        mylist <- strsplit(readLines(fc), "\t")
        close(fc)
        
        entry <- mylist[[grep(motif, mylist)]]
        tfs <- entry[2:length(entry)]
        plot_genes(tfs, F, paste0("tfs-of-", motif))
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
