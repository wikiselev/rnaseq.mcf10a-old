## RNA-Seq PIP3 signaling project

This is a collection of scripts that I used to analyze RNA-Seq time course data from MCF10A human breast cell lines upon EGF stimulation (postdoctoral project at the Babraham Institute, Cambridge, UK).

### External sources

1. [Draft of the manuscript](https://drive.google.com/folderview?id=0B9AEJU3ZybXIYkJ1T3JubFlOSWc&usp=sharing) (accessible only by collaborators).

2. Scripts for processing of the raw data can be found in `raw-processing` folder.

3. [Interactive RNA-Seq gene profiles](http://www.bioinformatics.babraham.ac.uk/shiny/kiselev-pip3-rna-seq-gene-profiles/)
(source files of this [Shiny app](http://shiny.rstudio.com/) can be found in `kiselev-pip3-rna-seq-gene-profiles` folder).

4. [ISMARA](http://ismara.unibas.ch/fcgi/mara) report can be viewed [here](http://lenoverelab.org/data/2015/kiselev/ismara_report_hg19/).

### main.R

To reproduce paper results one needs to run `main.R` which follows the paper chapters and where all the analysis steps are described. This file sources `functions.R` file which loads all required libraries and several function files (split by functionality):

* main-shortcuts.R
* read-process.R
* diff-expr.R
* ismara.R
* clustering.R
* plot.R
* go.R
* import-other-data.R
* alt-splicing.R
* linc-rna.R
* mirna.R

### extra-main.R

Some extra analysis, not included in the paper, can be performed by running `extra-main.R` which calls the following function files (split by functionality):

* alt-splicing.R
* splicing-on-cluster.R
* linc-rna.R
* mirna.R

### Data files

Data files will become available once the paper is published.
