#####
# SQL
#####

create_sqlite_database <- function() {
  file.remove("../rna-seq-media/data/rna-seq-db")
  sqldf("attach '../rna-seq-media/data/rna-seq-db' as new", drv = "SQLite")

  # gene_name_annotations table
  # ensembl_gene_id will be my key in the database
  d <- readRDS("../rna-seq-media/data/count_matrix_ann.rds")

  sqldf("create table gene_name_annotations as select * from d",
    dbname = "../rna-seq-media/data/rna-seq-db")

  d <- as.data.frame(readRDS("../rna-seq-media/data/diff_expr_wt.rds")$LRT)
  d$ensembl_gene_id <- rownames(d)
  sqldf("create table de_wt_time as select * from d",
    dbname = "../rna-seq-media/data/rna-seq-db")

  d <- as.data.frame(readRDS("../rna-seq-media/data/diff_expr_wt_ko.rds")$ko_vs_wt)
  d$ensembl_gene_id <- rownames(d)
  sqldf("create table de_wt_a66_time_cond as select * from d",
    dbname = "../rna-seq-media/data/rna-seq-db")

  d <- as.data.frame(readRDS("../rna-seq-media/data/diff_expr_0.rds")$A66)
  d$ensembl_gene_id <- rownames(d)
  sqldf("create table de_wt_a66_0_cond as select * from d",
    dbname = "../rna-seq-media/data/rna-seq-db")

  d <- as.data.frame(readRDS("../rna-seq-media/data/diff_expr_0.rds")$KI)
  d$ensembl_gene_id <- rownames(d)
  sqldf("create table de_wt_ki_0_cond as select * from d",
    dbname = "../rna-seq-media/data/rna-seq-db")

  d <- as.data.frame(readRDS("../rna-seq-media/data/diff_expr_0.rds")$PTEN)
  d$ensembl_gene_id <- rownames(d)
  sqldf("create table de_wt_pten_0_cond as select * from d",
    dbname = "../rna-seq-media/data/rna-seq-db")

  d <- as.data.frame(readRDS("../rna-seq-media/data/ismara_targets.rds"))
  sqldf("create table ismara_targets as select * from d",
    dbname = "../rna-seq-media/data/rna-seq-db")

  sqldf("select * from sqlite_master", dbname = "rna-seq-db")
}

get_gene_sets_from_database <- function() {
  sqldf(, dbname = "../rna-seq-media/data/rna-seq-db", drv = "SQLite")

  egf <- 
  c(sqldf(
    paste0("
      select ensembl_gene_id
      from de_wt_time
      where padj < ", padj
    )
  ))

  a66 <- 
  c(sqldf(
    paste0("
      select ensembl_gene_id
      from de_wt_a66_time_cond
      where padj < ", padj
    )
  ))

  ki <- 
  c(sqldf(
    paste0("
      select ensembl_gene_id
      from de_wt_ki_0_cond
      where padj < ", padj
    )
  ))

  pten <- 
  c(sqldf(
    paste0("
      select ensembl_gene_id
      from de_wt_pten_0_cond
      where padj < ", padj
    )
  ))


  a66.egf <-
  c(sqldf(
    paste0("
      select ensembl_gene_id
      from de_wt_a66_time_cond
      where padj < ", padj, " and ensembl_gene_id
      in (
        select ensembl_gene_id
        from de_wt_time
        where padj < ", padj, ")"
    )
  ))

  a66.negf <-
  c(sqldf(
    paste0("
      select ensembl_gene_id
      from de_wt_a66_time_cond
      where padj < ", padj, " and ensembl_gene_id
      not in (
        select ensembl_gene_id
        from de_wt_time
        where padj < ", padj, ")"
    )
  ))

  ki.egf <-
  c(sqldf(
    paste0("
      select ensembl_gene_id
      from de_wt_ki_0_cond
      where padj < ", padj, " and ensembl_gene_id
      in (
        select ensembl_gene_id
        from de_wt_time
        where padj < ", padj, ")"
    )
  ))

  ki.negf <-
  c(sqldf(
    paste0("
      select ensembl_gene_id
      from de_wt_ki_0_cond
      where padj < ", padj, " and ensembl_gene_id
      not in (
        select ensembl_gene_id
        from de_wt_time
        where padj < ", padj, ")"
    )
  ))

  pten.egf <-
  c(sqldf(
    paste0("
      select ensembl_gene_id
      from de_wt_pten_0_cond
      where padj < ", padj, " and ensembl_gene_id
      in (
        select ensembl_gene_id
        from de_wt_time
        where padj < ", padj, ")"
    )
  ))

  pten.negf <-
  c(sqldf(
    paste0("
      select ensembl_gene_id
      from de_wt_pten_0_cond
      where padj < ", padj, " and ensembl_gene_id
      not in (
        select ensembl_gene_id
        from de_wt_time
        where padj < ", padj, ")"
    )
  ))

  a66.ki <-
  c(sqldf(
    paste0("
      select ensembl_gene_id
      from de_wt_a66_time_cond
      where padj < ", padj, " and ensembl_gene_id
      in (
        select ensembl_gene_id
        from de_wt_ki_0_cond
        where padj < ", padj, ")"
    )
  ))

  a66.pten <-
  c(sqldf(
    paste0("
      select ensembl_gene_id
      from de_wt_a66_time_cond
      where padj < ", padj, " and ensembl_gene_id
      in (
        select ensembl_gene_id
        from de_wt_pten_0_cond
        where padj < ", padj, ")"
    )
  ))

  ki.pten <-
  c(sqldf(
    paste0("
      select ensembl_gene_id
      from de_wt_ki_0_cond
      where padj < ", padj, " and ensembl_gene_id
      in (
        select ensembl_gene_id
        from de_wt_pten_0_cond
        where padj < ", padj, ")"
    )
  ))

  ki.pten.egf <-
  c(sqldf("select ensembl_gene_id
    from de_wt_ki_0_cond
    where padj < 0.01 and ensembl_gene_id
    in (
      select ensembl_gene_id
      from de_wt_pten_0_cond
      where padj < 0.01
    ) and ensembl_gene_id
    in (
      select ensembl_gene_id
      from de_wt_time
      where padj < 0.01
    )"))

  ki.pten.negf <-
  c(sqldf("select ensembl_gene_id
    from de_wt_ki_0_cond
    where padj < 0.01 and ensembl_gene_id
    in (
      select ensembl_gene_id
      from de_wt_pten_0_cond
      where padj < 0.01
    ) and ensembl_gene_id
    not in (
      select ensembl_gene_id
      from de_wt_time
      where padj < 0.01
    )"))
}

#####
# SQL
#####
