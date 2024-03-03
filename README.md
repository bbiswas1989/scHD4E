Install the scHD4E package using the following: devtools::install_github("bbiswas1989/scHD4E")
Installation of scHD4E depends on the installation of some dependent package: scDEA, limma, ROSeq, Seurat, MAST, data.table, plyr, dplyr, gmodels, SwarnSeq.
These packages also depends on many other packages. All the dependencies should be install correctly before installing scHD4E. 
We suggests the use of R software with version 4.1.1 for the implementation of scHD4E.

scDEA (https://github.com/Zhangxf-ccnu/scDEA): The dependencies packages of scDEA are BPSC, DEsingle, DESeq2, edgeR, MAST, monocle, scDD, limma, Seuart, zingeR, scater,
SingleCellExperiment and aggregation. The guidelines of these package installation includs in scDEA in github.

ROSeq (https://github.com/krishan57gupta/ROSeq): Installation guides are here.

SwarnSeq (https://rdrr.io/github/sam-uofl/SwarnSeq/): Installation code are included in the link.

Package devtools are required to install many dependent package.

Some package installation code are given below:

Install BiocManager:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")

Install edgeR:
BiocManager::install("edgeR")
library(edgeR)

Install Limma:
BiocManager::install("limma")
library(limma)

Install DESeq2:
BiocManager::install("DESeq2")
library(DESeq2)

Install scDD:
BiocManager::install("scDD")
library(scDD)

Install BPSC:
devtools::install_github("nghiavtr/BPSC")
library(BPSC)

BiocManager::install("DEsingle")
library(DEsingle)

BiocManager::install("MAST")
library(MAST)

BiocManager::install("monocle")
library(monocle)

BiocManager::install("Seurat")
library(Seurat)

install.packages("remotes")
remotes::install_github("sam-uofl/SwarnSeq")

install.packages("devtools")
library("devtools")

devtools::install_github("statOmics/zingeR")
library(zingeR)

BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)

BiocManager::install("scater")
library(scater)

Dependencies must be installed correctly before the installation of scHD4E package, otherwise it may does not work properly.

install.packages("aggregation")
library(aggregation)

devtools::install_github("Zhangxf-ccnu/scDEA")
library("scDEA")

install_github('krishan57gupta/ROSeq')
library(ROSeq)
