Install the scHD4E package using the following: devtools::install_github("bbiswas1989/scHD4E")
Installation of scHD4E depends on the installation of some dependent packages: scDEA, limma, ROSeq, Seurat, MAST, data.table, plyr, dplyr, gmodels, SwarnSeq.
These packages also depend on many other packages. All the dependencies should be installed correctly before installing scHD4E. 
We have used the R software with version 4.1.1 for the development of the scHD4E package.

scDEA (https://github.com/Zhangxf-ccnu/scDEA): The dependencies packages of scDEA are BPSC, DEsingle, DESeq2, edgeR, MAST, monocle, scDD, limma, Seuart, zingeR, scater,
SingleCellExperiment and aggregation. The guidelines of these package installations include in scDEA in github.

ROSeq (https://github.com/krishan57gupta/ROSeq): Installation guides are here.

SwarnSeq (https://rdrr.io/github/sam-uofl/SwarnSeq/): Installation code are included in the link.

Package devtools are required to install many dependent package.

Some package installation codes are given below:

Install BiocManager:

if (!require("BiocManager", quietly = TRUE))

    install.packages("BiocManager")
    
BiocManager::install(version = "3.14")

BiocManager::install("edgeR")

BiocManager::install("limma")

BiocManager::install("DESeq2")

BiocManager::install("scDD")

devtools::install_github("nghiavtr/BPSC")

BiocManager::install("DEsingle")

BiocManager::install("MAST")

BiocManager::install("monocle")

BiocManager::install("Seurat")

install.packages("remotes")

remotes::install_github("sam-uofl/SwarnSeq")

install.packages("devtools")

devtools::install_github("statOmics/zingeR")

BiocManager::install("SingleCellExperiment")

BiocManager::install("scater")

install.packages("aggregation")

devtools::install_github("Zhangxf-ccnu/scDEA")

install_github('krishan57gupta/ROSeq')

Dependencies must be installed correctly before the installation of scHD4E package, otherwise it may does not work properly.


Locally run the shiny application by downloading the R programm file from 
the following website link: https://github.com/bbiswas1989/scHD4E-Shiny/blob/main/R/scHD4E-Shiny.R
