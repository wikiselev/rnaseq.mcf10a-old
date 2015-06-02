diff_expr_pairwise <- function(cond1, cond2, time1, time2) {
    files <- list.files("../pip3-rna-seq-input/htseq-read-count-deseq/")
    files <- files[grepl(paste0(cond1, "_", ".*", "_", time1, "_"), files) | grepl(paste0(cond2, "_", ".*", 
        "_", time2, "_"), files)]
    condition <- sapply(strsplit(files, "\\_"), "[[", 1)
    time <- sapply(strsplit(files, "\\_"), "[[", 3)
    condition.levels = c(paste(cond1, time1, sep = "_"), paste(cond2, time2, sep = "_"))
    samples <- data.frame(sample.name = files, file.name = files, condition = paste(condition, time, sep = "_"))
    
    # change levels
    samples$condition <- factor(samples$condition, levels = condition.levels)
    
    # create a DESeqDataSet with interaction term in the design formula
    cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples, directory = "../pip3-rna-seq-input/htseq-read-count-deseq/", 
        design = ~condition)
    
    # perform the main computation
    cds <- DESeq(cds)
    
    saveRDS(results(cds), paste0("../pip3-rna-seq-output/rds/diff-expr-", cond1, "-", time1, "-", cond2, "-", 
        time2, ".rds"))
}

diff_expr_pairwise_mirna <- function(cond1, cond2, time1, time2) {
    matrix <- readRDS("../pip3-rna-seq-output/rds/count-matrix-mirna.rds")
    matrix <- matrix[, grepl(paste0(cond1, "_", ".*", "_", time1), colnames(matrix)) | grepl(paste0(cond2, 
        "_", ".*", "_", time2), colnames(matrix))]
    
    condition <- sapply(strsplit(colnames(matrix), "\\_"), "[[", 1)
    condition.levels = c(paste(cond1, time1, sep = "_"), paste(cond2, time2, sep = "_"))
    
    time <- sapply(strsplit(colnames(matrix), "\\_"), "[[", 3)
    # time.levels <- c(0, 15, 40, 90, 180, 300)
    colData <- data.frame(sample.name = colnames(matrix), condition = paste(condition, time, sep = "_"))
    # change levels
    colData$condition <- factor(colData$condition, levels = condition.levels)
    
    # create a DESeqDataSet with interaction term in the design formula
    dds <- DESeqDataSetFromMatrix(countData = matrix, colData = colData, design = ~condition)
    
    # perform the main computation
    dds <- DESeq(dds)
    
    saveRDS(results(dds), paste0("../pip3-rna-seq-output/rds/diff-expr-mirna-", cond1, "-", time1, "-", cond2, 
        "-", time2, ".rds"))
}

diff_expr_time_course_mirna <- function(cond) {
    matrix <- readRDS("../pip3-rna-seq-output/rds/count-matrix-mirna.rds")
    matrix <- matrix[, grepl(cond, colnames(matrix))]
    
    time <- sapply(strsplit(colnames(matrix), "\\_"), "[[", 3)
    time.levels <- c(0, 15, 40, 90, 180, 300)
    colData <- data.frame(sample.name = colnames(matrix), time = time)
    # change levels
    colData$time <- factor(colData$time, levels = time.levels)
    
    # create a DESeqDataSet with interaction term in the design formula
    dds <- DESeqDataSetFromMatrix(countData = matrix, colData = colData, design = ~time)
    
    # perform the main computation
    dds <- DESeq(dds)
    
    ddsLRT <- nbinomLRT(dds, reduced = ~1)
    saveRDS(results(ddsLRT), paste0("../pip3-rna-seq-output/rds/diff-expr-mirna-", cond, ".rds"))
}

diff_expr_time_course <- function(cond) {
    files <- list.files("../pip3-rna-seq-input/htseq-read-count-deseq/")
    files <- files[grepl(paste0(cond, "_"), files)]
    time <- sapply(strsplit(files, "\\_"), "[[", 3)
    time.levels <- c("0", "15", "40", "90", "180", "300")
    samples <- data.frame(sample.name = files, file.name = files, time = time)
    # change levels
    samples$time <- factor(samples$time, levels = time.levels)
    
    # create a DESeqDataSet with interaction term in the design formula
    cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples, directory = "../pip3-rna-seq-input/htseq-read-count-deseq/", 
        design = ~time)
    
    # perform the main computation
    cds <- DESeq(cds)
    
    # using a new feature in DESeq2 package use maximum likelihood test to find diff. expressed genes in WT
    # time course
    cdsLRT <- nbinomLRT(cds, reduced = ~1)
    saveRDS(results(cdsLRT), paste0("../pip3-rna-seq-output/rds/diff-expr-", cond, ".rds"))
}

diff_expr_two_time_courses <- function(cond1, cond2) {
    files <- list.files("../pip3-rna-seq-input/htseq-read-count-deseq/")
    files <- files[grepl(paste(c(paste0(cond1, "_"), paste0(cond2, "_")), collapse = "|"), files)]
    time <- sapply(strsplit(files, "\\_"), "[[", 3)
    condition <- sapply(strsplit(files, "\\_"), "[[", 1)
    samples <- data.frame(sample.name = files, file.name = files, condition = condition, time = time)
    
    # change levels
    samples$condition <- factor(samples$condition, levels = c(cond1, cond2))
    samples$time <- factor(samples$time, levels = c("0", "15", "40", "90", "180", "300"))
    
    # create a DESeqDataSet with interaction term in the design formula
    cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples, directory = "../pip3-rna-seq-input/htseq-read-count-deseq/", 
        design = ~condition + time + condition:time)
    
    # perform the main computation
    cds <- DESeq(cds)
    resultsNames(cds)
    
    # likelihood ratio test
    cdsLRT <- nbinomLRT(cds, reduced = ~condition + time)
    saveRDS(results(cdsLRT), paste0("../pip3-rna-seq-output/rds/diff-expr-", cond1, "-", cond2, ".rds"))
}

diff_expr_two_time_courses_cond <- function(cond1, cond2) {
    files <- list.files("../pip3-rna-seq-input/htseq-read-count-deseq/")
    files <- files[grepl(paste(c(paste0(cond1, "_"), paste0(cond2, "_")), collapse = "|"), files)]
    time <- sapply(strsplit(files, "\\_"), "[[", 3)
    condition <- sapply(strsplit(files, "\\_"), "[[", 1)
    samples <- data.frame(sample.name = files, file.name = files, condition = condition, time = time)
    
    # change levels
    samples$condition <- factor(samples$condition, levels = c(cond1, cond2))
    samples$time <- factor(samples$time, levels = c("0", "15", "40", "90", "180", "300"))
    
    # create a DESeqDataSet with interaction term in the design formula
    cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples, directory = "../pip3-rna-seq-input/htseq-read-count-deseq/", 
        design = ~condition + time + condition:time)
    
    # perform the main computation
    cds <- DESeq(cds)
    resultsNames(cds)
    
    # likelihood ratio test
    cdsLRT <- nbinomLRT(cds, reduced = ~time)
    saveRDS(results(cdsLRT), paste0("../pip3-rna-seq-output/rds/diff-expr-", cond1, "-", cond2, "-cond.rds"))
}

get_diff_expr <- function(name, padj) {
    d <- readRDS(paste0("../pip3-rna-seq-output/rds/diff-expr-", name, ".rds"))
    return(d[!is.na(d$padj) & d$padj < padj, ])
}


get_diff_expr_time_course_pairwise <- function(cond, padj) {
    t_15 <- rownames(get_diff_expr(paste(cond, "0", cond, "15", sep = "-"), padj))
    t_40 <- rownames(get_diff_expr(paste(cond, "0", cond, "40", sep = "-"), padj))
    t_90 <- rownames(get_diff_expr(paste(cond, "0", cond, "90", sep = "-"), padj))
    t_180 <- rownames(get_diff_expr(paste(cond, "0", cond, "180", sep = "-"), padj))
    t_300 <- rownames(get_diff_expr(paste(cond, "0", cond, "300", sep = "-"), padj))
    
    return(unique(c(t_15, t_40, t_90, t_180, t_300)))
}

venn <- function(list, is.weights, name) {
    v <- Venn(list)
    pdf(file = paste0("../pip3-rna-seq-output/figures/venn-", name, ".pdf"))
    plot(v, doWeights = is.weights)
    dev.off()
}

venn_ellipses <- function(list, is.weights, name) {
    v <- Venn(list)
    pdf(file = paste0("../pip3-rna-seq-output/figures/venn-", name, ".pdf"))
    plot(v, doWeights = is.weights, type = "ellipses")
    dev.off()
} 
