################
# ANOVA analysis
################

two_way_anova <- function(dat, ret.val) {
  res <- anova(lm(value ~ cond + time + cond * time, dat))
  if(ret.val == "int") {
    return(res[[5]][3])
  }
  if(ret.val == "cond") {
    return(res[[5]][1])
  }
  if(ret.val == "time") {
    return(res[[5]][2])
  }
}

anova_results <- function(cond1, cond2, genes) {
  anova_res_int <- count.data.table[gene_id %in% genes & (cond == cond1 | cond == cond2),
    list(pval =
      two_way_anova(data.frame(value = value, cond = cond, time = time), "int")),
    by = gene_id]
  anova_res_cond <- count.data.table[gene_id %in% genes & (cond == cond1 | cond == cond2),
    list(pval =
      two_way_anova(data.frame(value = value, cond = cond, time = time), "cond")),
    by = gene_id]
  # anova_res_time <- count.data.table[gene_id %in% genes & (cond == cond1 | cond == cond2),
  #   list(pval =
  #     two_way_anova(data.frame(value = value, cond = cond, time = time), "time")),
  #   by = gene_id]

  anova_res_int$padj <- p.adjust(anova_res_int$pval, method = "BH")
  anova_res_cond$padj <- p.adjust(anova_res_cond$pval, method = "BH")
  # anova_res_time$padj <- p.adjust(anova_res_time$pval, method = "BH")
  return(list(int = anova_res_int, cond = anova_res_cond))#, time = anova_res_time))
}

nova <- function(genes) {
  res.a66 <- anova_results("wt", "a66", genes)
  res.pten <- anova_results("wt", "pten", genes)
  res.ki <- anova_results("wt", "ki", genes)

  res <- data.table(gene_id = res.a66$int[,gene_id])
  res$a66_int <- res.a66$int[,padj]
  res$a66_cond <- res.a66$cond[,padj]

  res$pten_int <- res.pten$int[,padj]
  res$pten_cond <- res.pten$cond[,padj]

  res$ki_int <- res.ki$int[,padj]
  res$ki_cond <- res.ki$cond[,padj]
  return(res)
}

get_sig_genes <- function(res) {
  # significant interaction effect between wt and a66
  a66_genes1 <- res[a66_int <= 1e-2, gene_id]
  # insignificant interaction between wt and a66 but significant condition effect
  a66_genes2 <- res[a66_cond <= 1e-3 & a66_int > 1e-2, gene_id]
  # significant interaction effect between wt and pten
  pten_genes1 <- res[pten_int <= 1e-2, gene_id]
  # insignificant interaction between wt and pten but significant condition effect
  pten_genes2 <- res[pten_cond <= 1e-3 & pten_int > 1e-2, gene_id]
  # significant interaction effect between wt and ki
  ki_genes1 <- res[ki_int <= 1e-2, gene_id]
  # insignificant interaction between wt and ki but significant condition effect
  ki_genes2 <- res[ki_cond <= 1e-3 & ki_int > 1e-2, gene_id]
  return(unique(c(a66_genes1, a66_genes2, pten_genes1, pten_genes2, ki_genes1, ki_genes2)))
}

################
# ANOVA analysis
################
