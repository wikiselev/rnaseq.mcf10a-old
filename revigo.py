#!/usr/bin/python
 
"""
- Submit example data to REVIGO server (http://revigo.irb.hr/)
- Download and run R script for creating the treemap
- Download and run R script for creating the scatterplot
 
Creates files:
    treemap.R, treemap.Rout, revigo_treemap.pdf
    scatter.R, scatter.Rout, revigo_scatter.pdf
"""
 
import os, sys
import urllib
import mechanize
import shutil

args = sys.argv

os.chdir("../pip3-rna-seq-output/GO/" + args[2])
 
url = "http://revigo.irb.hr/"
 
# RobustFactory because REVIGO forms not well-formatted
br = mechanize.Browser(factory=mechanize.RobustFactory())
 
# For actual data, use open('mydata.txt').read()
# br.open(os.path.join(url, 'examples', 'example1.txt'))
# txt = br.response().read()
txt = open(args[1]).read()

# directory = args[1].split(".")[0]
# directory = "../revigo/" + args[2] + "/" + directory
# # shutil.rmtree(directory)
# os.makedirs(directory)

# Encode and request
data = {'inputGoList': txt}
br.open(url, data=urllib.urlencode(data))

# Choose HomoSapiens in the drop down menu
br.select_form(nr=0)
control = br.form.find_control("goSizes")
# loop through items to find the match
for item in control.items:
  if item.name == "9606":

    # it matches, so select it
    item.selected = True
 
# Submit form
br.select_form(name="submitToRevigo")
response = br.submit()
 
# Exact string match on the url for getting the R treemap script
br.follow_link(url="export.jsp?table=1")
with open('REVIGO.csv', 'w') as f:
    f.write(br.response().read())

br.back()

br.follow_link(url="toR_treemap.jsp?table=1")
with open('treemap.R', 'w') as f:
    f.write(br.response().read())
 
# Create PDFs
# os.chdir(args[2])
os.system('Rscript treemap.R')
