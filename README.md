RNA-Seq PIP3 signaling project
=============

This is a collection of scripts that I used to analyze RNA-Seq time course data from MF10a human breast cell lines upon EGF stimulation (postdoctoral project at the Babraham Institute, Cambridge, UK).

The paper is in preparation at the moment. Draft is available [here](https://drive.google.com/folderview?id=0B9AEJU3ZybXIYkJ1T3JubFlOSWc&usp=sharing) (accessible only by collaborators).

Raw ISMARA report can be viewed [here](http://lenoverelab.org/data/2015/kiselev/ismara_report_hg19/).

Averaged ISMARA report can be viewed [here](http://lenoverelab.org/data/2015/kiselev/averaged_report_hg19/).

The main script is `main.R` which follows the paper chapters and where all the analysis steps are described. This file sources `functions.R` file which loads all required libraries and several function files (split by functionality):

main.R

* functions.R

  * read-process.R

  * diff-expr.R

  * const-mut.R

  * time-courses.R

  * ismara.R

  * sql.R

  * plot.R

  * svd.R

  * anova.R

  * go.R
