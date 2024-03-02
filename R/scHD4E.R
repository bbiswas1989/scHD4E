# This package has been developed to implement the scHD4E method.
#
# The main function of this package is scHD4E.
# It contains the object data matrix, condition (case/control) and ncluster (number of cluster)
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

###########################################
################scHD4E#####################
###########################################

###### Import these Packages########
library("scDEA")
library(MAST)
library(data.table)
library(plyr)
library(dplyr)
library(gmodels)
library(SwarnSeq)
library(ROSeq)
library(edgeR)
library(limma)

#####Main function for scHD4E method#######

scHD4E<-function(scData,condition,ncluster){

  #This looping is necessary when a number of genes of a dataset contains all zeros across cell#

  index<-c()
  for (i in 1:dim(scData)[1]){
    index[i]=length(which(as.numeric(scData[i,])!=0))/length(as.numeric(scData[i,]))
  }
  keep1<-which(index>=0.001)
  filterData<-scData[keep1,]

  ############ROSeq##########
  pred.ROSeq <- Sys.time()
  output<-ROSeq(countData=filterData, condition = condition, numCores=1)
  data_ROSeq<-data.frame(output)

  #### The program of the next section is essential when a method produce missing p-values##
  dat1<-data.frame(scData[,1:2])
  library(tibble)
  dat1 = rownames_to_column(dat1, "Plants")
  dat2 = rownames_to_column(data_ROSeq, "Plants")
  library(dplyr)
  dat = full_join(dat1, dat2, )
  dat = dat %>% replace(is.na(.), 0.90)   ####when a gene produce missing value we replace it by 0.90. Thats means we assume it is non-differential. This is nessary to avoid the disruption of algorithm#####
  #####################################

  ROSeq<-dat$pVals
  end.ROSeq<- Sys.time()
  time.ROSeq<- difftime(end.ROSeq, pred.ROSeq, units = "mins")
  cat("Run time for ROSeq: ", time.ROSeq, "min","\n")

  #########Seurat##############
  pred.Seurat <- Sys.time()
  Pvals_Seurat_full<-scDEA_individual_methods(
    raw.count=scData,
    cell.label=condition,
    is.normalized = FALSE,
    verbose = TRUE,
    BPSC = F,
    DEsingle = F,
    DESeq2 = F,
    edgeR = F,
    MAST = F,
    monocle = F,
    scDD = F,
    Ttest = F,
    Wilcoxon = F,
    limma = F,
    Seurat = TRUE,
    zingeR.edgeR = F,
    Seurat.normalize = "CPM",
    Seurat.method = "bimod")
  Seurat<-Pvals_Seurat_full
  end.Seurat<- Sys.time()
  time.Seurat <- difftime(end.Seurat, pred.Seurat, units = "mins")
  cat("Run time for Seurat: ", time.Seurat, "min","\n")

  ###########Limma-Voom##################
  pred.Limma <- Sys.time()
  design_full<- model.matrix(~condition)
  nf_full<- calcNormFactors(scData)
  v_full<- voom(scData, design_full, lib.size=colSums(scData)*nf_full,
                normalize.method="quantile", plot=F)
  fit.voom_full<- lmFit(v_full, design_full)
  fit.voom_full<- eBayes(fit.voom_full)
  Limma<-fit.voom_full$p.value[,2]
  end.Limma<- Sys.time()
  time.Limma <- difftime(end.Limma, pred.Limma, units = "mins")
  cat("Run time for Limma: ", time.Limma, "min","\n")

  ##############TPMM##############
  pred.TPMM <- Sys.time()
  fData_full<- data.frame(primerid=rownames(scData))

  ##### next line code is used to determine the number of cluster using SwarnSeq R package##
  #results <- optimcluster(CountData=as.matrix(scDataB), n = 10, seed = 108, Threshold = 0.3, plot = T)
  set.seed(30)

  #####We include clustering effect instead of batch effect due to the absence of batch effect for most of the data#####

  km.res_full<- kmeans(t(scData), ncluster, nstart = 10)
  PoolID_full<-as.vector(km.res_full$cluster)
  cData_full<- data.frame(wellKey=colnames(scData),PoolID_full)
  count_full<-as.matrix(scData)
  scaRaw_full<- FromMatrix(count_full, cData_full, fData_full,check_sanity =F)

  lmer.output_full<- zlm(~condition+(1|PoolID_full), scaRaw_full, method='bayesglm', ebayes=T)
  lmer.lr_full<- lrTest(lmer.output_full, 'condition')
  mixed_full= cbind(lmer.lr_full[,,'Pr(>Chisq)'], lmer.lr_full[,,'lambda'])
  colnames(mixed_full) = c('gaussian_mixed_pvalue','binomial_mixed_pvalue','combine_mixed_pvalue','gaussian_mixed_lambda','binomial_mixed_lambda','combine_mixed_lambda')

  data_TPMM_full<-data.frame(mixed_full)
  TPMM<-data_TPMM_full$combine_mixed_pvalue

  end.TPMM<- Sys.time()
  time.TPMM <- difftime(end.TPMM, pred.TPMM, units = "mins")
  cat("Run time for TPMM: ", time.TPMM, "min","\n")


  #######scDEA##############
  pred.scDEA <- Sys.time()
  Pvals_scDEA<-scDEA_individual_methods(
    raw.count=scData,
    cell.label=condition,
    is.normalized = FALSE,
    verbose = TRUE,
    BPSC = TRUE,
    DEsingle = TRUE,
    DESeq2 = TRUE,
    edgeR = TRUE,
    MAST = TRUE,
    monocle = TRUE,
    scDD = TRUE,
    Ttest = TRUE,
    Wilcoxon = TRUE,
    limma = TRUE,
    Seurat = TRUE,
    zingeR.edgeR = TRUE,
    BPSC.coef = 2,
    BPSC.normalize = "CPM",
    BPSC.parallel = TRUE,
    DEsingle.parallel = TRUE,
    DEsingle.normalize = "CPM",
    DESeq2.test = "LRT",
    DESeq2.parallel = TRUE,
    DESeq2.beta.prior = TRUE,
    DESeq2.fitType = "parametric",
    DESeq2.normalize = "CPM",
    edgeR.Test = "QLFT",
    edgeR.normalize = "TMM",
    limma.method.fit = "ls",
    limma.trend = TRUE,
    limma.robust = TRUE,
    limma.normalize = "CPM",
    Seurat.normalize = "CPM",
    Seurat.method = "bimod",
    MAST.method = "bayesglm",
    MAST.normalize = "CPM",
    MAST.parallel = TRUE,
    monocle.cores = 1,
    monocle.normalize = "CPM",
    scDD.alpha1 = 0.01,
    scDD.mu0 = 0,
    scDD.s0 = 0.01,
    scDD.a0 = 0.01,
    scDD.b0 = 0.01,
    scDD.normalize = "CPM",
    scDD.permutation = 0,
    Ttest.normalize = "CPM",
    Wilcoxon.normalize = "CPM",
    zingeR.edgeR.normalize = "CPM",
    zingeR.edgeR.maxit.EM = 100
  )

  combination.Pvals<- lancaster.combination(Pvals_scDEA, weight = TRUE, trimmed = 0.2)
  adjusted.Pvals<- scDEA.p.adjust(combination.Pvals, adjusted.method = "BH")
  scDEA<-data.frame(Pvals_scDEA,combination.Pvals,adjusted.Pvals)

  end.scDEA<- Sys.time()
  time.scDEA <- difftime(end.scDEA, pred.scDEA, units = "mins")
  cat("Run time for scDEA: ", time.scDEA, "min","\n")

  #############scHD4E############
  pred.HD4E<- Sys.time()
  Pvals_HN_full<-as.matrix(data.frame(Seurat,Limma,TPMM,ROSeq))
  Pvals_HN_full= Pvals_HN_full%>% replace(is.na(.), 0.20)      # We impute the missing p-values by 0.20 if produce any missing value(rare). This indicates that the genes are non-differentially expressed##

  combination.Pvals_HN_full<- lancaster.combination(Pvals_HN_full, weight = TRUE, trimmed = 0.2)
  adjusted.Pvals_HN_full<- scDEA.p.adjust(combination.Pvals_HN_full, adjusted.method = "BH")
  end.HD4E<- Sys.time()
  time.HD4E1<- difftime(end.HD4E, pred.HD4E, units = "mins")
  time.HD4E<-time.Seurat+time.Limma+time.TPMM+time.ROSeq+time.HD4E1
  cat("Run time for HD4E: ", time.HD4E, "min","\n")

  data_HN_full<-data.frame(Pvals_HN_full,combination.Pvals_HN_full,adjusted.Pvals_HN_full,scDEA) ####Data frame of p-values including our proposed method#############
  return(data_HN_full)
}

######Curated scRNA-seq data######

#scDataB<-read.csv("F:\\PhD Folder\\Analysis for HD4E for GSE77288\\GSE81608_Xin_Top3000.csv",header=T,row.names=1)

#######Group########
#condition<- factor(c(rep("1",229),rep("2",161)))

#scData<-scDataB
#Pvals_full<-scHD4E(scData,condition,ncluster=5)
#ncluster=5 for GSE81608 dataset.
#ncluster is selected with the optimcluster function implementation under SwarnSeq package
#results <- optimcluster(CountData=as.matrix(scDataB), n = 10, seed = 108, Threshold = 0.3, plot = T)
