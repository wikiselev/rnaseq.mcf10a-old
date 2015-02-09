rm -r ../rna-seq-media/figures-new-png/
mkdir ../rna-seq-media/figures-new-png/
sips -s format png ../rna-seq-media/figures_new/*.pdf --out ../rna-seq-media/figures-new-png/

# Rscript revigo.R