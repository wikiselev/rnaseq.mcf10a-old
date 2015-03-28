miRTarBase <- function() {
        d <- read.csv("../pip3-rna-seq-input/miRNA/hsa_MTI.csv")
        saveRDS(as.data.table(d), "../pip3-rna-seq-output/rds/miRTarBase.rds")
}

get_miRTarBase <- function() {
        return(readRDS("../pip3-rna-seq-output/rds/miRTarBase.rds"))
}

mirna_read_process <- function() {
        # Read in the raw data
        count.matrix <- read.delim("../pip3-rna-seq-input/miRNA/summarised_mirna_counts.txt",
                                   row.names = 1)
        ann <- read.table("../pip3-rna-seq-input/miRNA/mirna_ann.txt")
        
        colnames(count.matrix) <- ann[,2][match(colnames(count.matrix), ann[,1])]
        count.matrix <- count.matrix[,!is.na(colnames(count.matrix))]
        
        # Remove miRNA which are not expressed in all conditions
        count.matrix <- as.matrix(count.matrix[rowSums(count.matrix) > 0, ])
        
        saveRDS(count.matrix, "../pip3-rna-seq-output/rds/count-matrix-mirna.rds")
        
        # Normalize by total read number
        cols <- colnames(count.matrix)
        count.matrix <- count_matrix_norm_matrix(count.matrix)
        colnames(count.matrix) <- cols
        
        saveRDS(count.matrix, "../pip3-rna-seq-output/rds/count-matrix-mirna-norm.rds")
        
        cor_heatmap(count.matrix, "cor-heatmap-mirna.pdf")
        
        conditions <- c("wt", "ko")
        
        res <- average_count_matrix(count.matrix, conditions)
        
        saveRDS(res[[1]], "../pip3-rna-seq-output/rds/count-matrix-mirna-av.rds")
        saveRDS(res[[2]], "../pip3-rna-seq-output/rds/count-matrix-mirna-sd.rds")
        
        plot.data.all <- create_plot_table(res[[1]], res[[2]])
        
        saveRDS(plot.data.all, "../pip3-rna-seq-output/rds/plot-time-courses-all-mirna.rds")
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

