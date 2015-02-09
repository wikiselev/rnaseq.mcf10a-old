get_constitutive_mutations <- function(padj) {
  d <- get_diff_expr("wt-0-ki-0", padj)
  ki <- data.frame(id = rownames(d), log2F = d$log2FoldChange, padj = d$padj)

  d <- get_diff_expr("wt-0-pten-0", padj)
  pten <- data.frame(id = rownames(d), log2F = d$log2FoldChange, padj = d$padj)

  d <- constitutive_mutations_data_table(ki, pten)
  return(d)
}

constitutive_mutations_data_table <- function(ki, pten) {
  d <- merge(ki, pten, by = "id", all = T)
  colnames(d) <- c("id", "log2F.ki", "padj.ki", "log2F.pten", "padj.pten")
  d <- as.data.table(d)
  d[,egf:=0]

  egf.genes <- rownames(get_diff_expr("wt", 0.01))

  d[id %in% egf.genes, egf:=1]
  d[is.na(log2F.ki), log2F.ki:=0]
  d[is.na(padj.ki), padj.ki:=1]
  d[is.na(log2F.pten), log2F.pten:=0]
  d[is.na(padj.pten), padj.pten:=1]
  return(d)
}

constitutive_mutations_sum_up <- function(d1, name) {
  d1[ki > 0, ki:=1]
  d1[ki < 0, ki:=-1]
  d1[pten > 0, pten:=1]
  d1[pten < 0, pten:=-1]
  
  t <- d1[,list(gene.num = length(id)),by=c("ki","pten","egf")]

  t$ki <- as.character(t$ki)
  t$pten <- as.character(t$pten)
  t$egf <- as.character(t$egf)

  t[ki == "1", ki:="up"]
  t[ki == "-1", ki:="down"]
  t[ki == "0", ki:="none"]  
  t[pten == "1", pten:="up"]
  t[pten == "-1", pten:="down"]
  t[pten == "0", pten:="none"]  
  t[egf == "1", egf:="yes"]
  t[egf == "0", egf:="no"]

  t[,ki_pten:=paste(t$ki, t$pten, sep = "_")]
  t$ki_pten = factor(t$ki_pten,
    levels = c("down_down", "up_up", "down_up", "up_down", "down_none",
    "up_none", "none_down", "none_up"))

  p <- ggplot(t, aes(ki_pten, gene.num)) +
        geom_bar(stat = "identity", position = "dodge", aes(fill = egf))

  pdf(file = paste0("../rna-seq-media/figures/constitutive_mutation_analysis_", name, ".pdf"),
    w = 10, h = 4)
  print(p)
  dev.off()

  return(t[,list(ki_pten, egf, gene.num)])
}

constitutive_mutations_go <- function(d, name) {
  GO(d[ki < 0 & pten < 0 & egf == 1, id], paste0("ki-pten-egf+", name))
  GO(d[ki < 0 & pten < 0 & egf == 0, id], paste0("ki-pten-egf-", name))
  GO(d[ki > 0 & pten > 0 & egf == 1, id], paste0("ki+pten+egf+", name))
  GO(d[ki > 0 & pten > 0 & egf == 0, id], paste0("ki+pten+egf-", name))
  GO(d[ki > 0 & pten < 0 & egf == 1, id], paste0("ki+pten-egf+", name))
  GO(d[ki > 0 & pten < 0 & egf == 0, id], paste0("ki+pten-egf-", name))
  GO(d[ki < 0 & pten > 0 & egf == 1, id], paste0("ki-pten+egf+", name))
  GO(d[ki < 0 & pten > 0 & egf == 0, id], paste0("ki-pten+egf-", name))

  GO(d[ki == 0 & pten > 0 & egf == 1, id], paste0("pten+egf+", name))
  GO(d[ki == 0 & pten > 0 & egf == 0, id], paste0("pten+egf-", name))
  GO(d[ki == 0 & pten < 0 & egf == 1, id], paste0("pten-egf+", name))
  GO(d[ki == 0 & pten < 0 & egf == 0, id], paste0("pten-egf-", name))

  GO(d[ki > 0 & pten == 0 & egf == 1, id], paste0("ki+egf+", name))
  GO(d[ki > 0 & pten == 0 & egf == 0, id], paste0("ki+egf-", name))
  GO(d[ki < 0 & pten == 0 & egf == 1, id], paste0("ki-egf+", name))
  GO(d[ki < 0 & pten == 0 & egf == 0, id], paste0("ki-egf-", name))
}

analyse_constitutive_mutations <- function(padj) {
  d <- get_constitutive_mutations(padj)
  constitutive_mutations_sum_up(d, padj)
  constitutive_mutations_go(d, padj)
}

process_const_mut_eff <- function(t, ki, pten, egf) {
  len <- dim(t)[1]
  t <- as.data.frame(t)
  t <- t[order(-t$log2F.eff),]
  if (len < 100) {
    end <- 0
  } else {
    end <- len%/%100 - 1
  }
  for (i in 0:end) {
    plot_genes_by_constitutive_eff(t[(i*100 + 1):((i + 1)*100),],
      paste0("ki", ki, "pten", pten, "egf", egf, i))
  }
  if (len > 100) {
    i <- i + 1
    plot_genes_by_constitutive_eff(t[(i*100+1):len,], 
      paste0("ki", ki, "pten", pten, "egf", egf, i))
  }
}

filtering_of_constitutive_effects <- function(padj) {
  d <- get_constitutive_mutations(padj)
  # filter out those genes where both mutations are effective
  d <- d[log2F.ki*log2F.pten != 0]
  # total effect of both mutations
  d[,log2F.eff:=abs(log2F.ki) + abs(log2F.pten)]

  # both mutations reduce gene expression at time 0
  process_const_mut_eff(d[log2F.ki < 0 & log2F.pten < 0 & egf==0], "-", "-", "-")
  process_const_mut_eff(d[log2F.ki < 0 & log2F.pten < 0 & egf==1], "-", "-", "+")
  # both mutations increase gene expression at time 0
  process_const_mut_eff(d[log2F.ki > 0 & log2F.pten > 0 & egf==0], "+", "+", "-")
  process_const_mut_eff(d[log2F.ki > 0 & log2F.pten > 0 & egf==1], "+", "+", "+")
}
