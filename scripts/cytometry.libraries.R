# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Libraries for cytometry analyses
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€ 

# Install packages
#
#source("https://bioconductor.org/biocLite.R")
#biocLite("flowClust")
#biocLite("openCyto")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("flowClust", "openCyto"))
# Load packages
library(flowViz)
library(flowClust)
library(flowMerge)
library(flowStats)
library(flowWorkspace)
library(openCyto)
library(flowCore)

#BiocInstaller::useDevel()
#biocLite("flowWorkspace")
#devtools::install_github("RGLab/ggcyto", ref = "trunk")

#library(ggcyto)