
alt_splicing <- function(cond1, cond2, time1, time2) {
  countFiles = list.files("../rna-seq-media/htseq_count_exons_new/",
    pattern="sam.txt$", full.names=TRUE)
  file.names = list.files("../rna-seq-media/htseq_count_exons_new/",
    pattern="sam.txt$", full.names=FALSE)
  countFiles <-
    countFiles[grepl(paste0(cond1, "_", ".*", "_", time1, "_"), countFiles) |
              grepl(paste0(cond2, "_", ".*", "_", time2, "_"), countFiles)]
  file.names <-
    file.names[grepl(paste0(cond1, "_", ".*", "_", time1, "_"), file.names) |
              grepl(paste0(cond2, "_", ".*", "_", time2, "_"), file.names)]
  samples <- paste(sapply(strsplit(file.names, "\\_"), "[[", 1),
    sapply(strsplit(file.names, "\\_"), "[[", 3),
    sapply(strsplit(file.names, "\\_"), "[[", 2), sep = "_")
  condition <- sapply(strsplit(file.names, "\\_"), "[[", 1)
  time <- sapply(strsplit(file.names, "\\_"), "[[", 3)

  sampleTable <- data.frame(row.names = samples,
    condition = condition, time = time)


  flattenedFile <- list.files("../rna-seq-media/htseq_count_exons_new/",
    pattern="gff$", full.names=TRUE)

  dxd <- DEXSeqDataSetFromHTSeq(
    countFiles,
    sampleData = sampleTable,
    design = ~ sample + exon + condition:exon,
    flattenedfile = flattenedFile
  )

  # head( counts(dxd), 5 )
  # split( seq_len(ncol(dxd)), colData(dxd)$exon )
  # head( featureCounts(dxd), 5 )
  # head( rowData(dxd), 3 )
  # sampleAnnotation( dxd )

  BPPARAM <- MulticoreParam(workers=4)
  dxd <- estimateSizeFactors( dxd )
  dxd <- estimateDispersions( dxd, BPPARAM=BPPARAM)
  dxd <- testForDEU( dxd, BPPARAM=BPPARAM)
  dxd <- estimateExonFoldChanges(dxd, BPPARAM=BPPARAM)

  dxr1 <- DEXSeqResults( dxd )

  saveRDS(dxr1, paste0("../rna-seq-media/data_new/exon_diff_expr_",
    cond1, "_", time1, "_", cond2, "_", time2, ".rds"))

  unlink(paste0("../rna-seq-media/DEXSeqReport_", cond2), recursive = T)
  DEXSeqHTML( dxr1, path = paste0("../rna-seq-media/DEXSeqReport_", cond2),
    FDR = 0.01, color=c("#FF000080", "#0000FF80") )
}

