# process mirna reads
mirna_read_process()
# import miRTarBase database with miRNA targets
miRTarBase()

# alternative splicing this should be done on a cluster because it takes several hours to finish please
# use script e.g. 'Rscript splicing_on_cluster.R pten' to compare wt and pten
alt_splicing("wt", "ki", 0, 0)
alt_splicing("wt", "pten", 0, 0)

# find out which of the genes are lincRNAs
linc_rna("mut")
linc_rna("a66")
# plot genes selected by Vero which might be regulated by lincRNAs (affected by mutations) need to rerun
# it on new lists - in the old ones some gene names corresponded to hundreds of ENSEMBL IDs!!!
plot_linc_rna() 
