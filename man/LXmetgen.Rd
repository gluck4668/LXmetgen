\name{LXmetgen}
\alias{LXmetgen}
\title{LXmetgen: gene and metabolite joint analysis}
\usage{
LXmetgen (gene_file,meta_file,group1,group2)
}
\description{
LXmetgen package can be used to analyze the joint-pathways between the genes and metabolites.
It cited several R packages publised by "https://www.metaboanalyst.ca/MetaboAnalyst".
}
\examples{

if(!requireNamespace("devtools"))
  install.packages("devtools")
library(devtools)

install_github("gluck4668/LXmetgen")
library(LXmetgen)

#---------------------------------
data(gene_data_example)
data(meta_data_example)

#------------------------------

setwd("D:/R-lin study/R packages/LXmetgen_example")

gene_file="gene_data.xlsx"
meta_file="meta_data.xlsx"

group1="Modle"
group2="Normal"

LXmetgen (gene_file,meta_file,group1,group2)


}
