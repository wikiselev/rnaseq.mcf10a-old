pdf_to_png <- function(files, name) {
	f <- files[grep(paste0("_", name, ".*.pdf"), files)]
	for(i in f) {
		new <- paste0("../pip3-rna-seq-output/revigo-png/", name, "/", strsplit(i, "\\/")[[1]][1], ".png")
		system(paste0("sips -s format png ../pip3-rna-seq-output/revigo/", i, " --out ", new))
	}
}

system("rm -r ../pip3-rna-seq-output/revigo-png/")
system("mkdir ../pip3-rna-seq-output/revigo-png/")
system("mkdir ../pip3-rna-seq-output/revigo-png/BP/")
system("mkdir ../pip3-rna-seq-output/revigo-png/MF/")
system("mkdir ../pip3-rna-seq-output/revigo-png/CC/")

files <- list.files("../pip3-rna-seq-output/revigo/", recursive = T)

pdf_to_png(files, "BP")
pdf_to_png(files, "MF")
pdf_to_png(files, "CC")
