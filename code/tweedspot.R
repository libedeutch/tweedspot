
library(data.table)
library(mgcv)
library(SpatialExperiment)
library(nnSVG)
library(tidyverse)
library(CTSV)
library(pscl)
library(DIRECT)
library(scran)
library(scuttle)

#load("~/OneDrive/cornell/simulData/gau.RData")

file_tot = list.files('/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/spot/zero_60/', ".RData")
#file_tot2 = list.files("~/OneDrive/cornell/simulMI/tweedspot/", ".RData")

tweedspot<- function(Y, W, spot.coor){
  
  #@`parameters`Y is nxG count matrix, W is all nx1, gamma_true is gamma, spot.coor is loc 
  
  #Y <- t(Y)
  loc <- spot.coor # coordinates are (1,1) (1,2) (1,3) .. 
  W = t(W)
  # make sure the sum of cell type proportions is equal to 1 in each spot.    
  W <- W / rowSums(W)
  # number of genes
  G <- ncol(Y)
  # number of spots
  n <- nrow(loc)
  # number of cell types
  K <- ncol(W)
  # normalize cell-type proportion matrix W to ensure the summation across cell types in one spot is equal to one.
  #W <- W / rowSums(W)
  # Center and normalize coordinates of spots to have mean zero and standard deviation one.
  S <- t(loc) - colMeans(loc)
  S <- t(S / apply(S, 1, sd))

  # DESIGN MATRIX 
  h1 <- S[,1]
  h2 <- S[,2]
  Tmp <- cbind(W * h1, W * h2, W)
  XY <- c("x", "y")
  Celltype <- paste('Celltype', 1:K, sep = '_')
  Celltype_XY<-do.call(paste0, expand.grid(Celltype, XY))
  colnames(Tmp)<-c(Celltype_XY, Celltype)
  #ell <- rowSums(Y) / median(rowSums(Y))
  # replace size factor 
  sce <- SingleCellExperiment(assays = list(counts = t(Y))) 
  sce <- scran::computeSumFactors(sce)
  ell <- sizeFactors(sce)
  K <- ncol(Tmp)/3
  
  #########################
  # TP GENES PER CELLTYPE #
  #########################
  
  
  pval = matrix(nrow = K, ncol = G)
  # Test 
  for (i in 1:ncol(Y)){
    #print(i)
    z <- Y[,i] 
    dat0<-cbind.data.frame(z, Tmp)
    #dat0<-cbind.data.frame(z, S)
    #celltype<- as.factor(max.col(W))
    #dat0$celltype<-celltype
    #fit0 <- gam(z ~ offset(log(ell)) + s(x, y, by = celltype), data = dat0, family ='tw')
    fit0 <- gam(z ~ - 1 + offset(log(ell)) + 
                  s(Celltype_1x, Celltype_1y) + 
                  s(Celltype_2x, Celltype_2y) + 
                  s(Celltype_3x, Celltype_3y) + 
                  Celltype_1 + Celltype_2 + Celltype_3  , 
                data = dat0, family ='tw', method = 'REML' )
    cc = summary(fit0)
    #pval[1,i] = cc$p.pv # or last column in cc$p.table
    pval[1:K,i] = cc$s.table[,'p-value']
    
   # save(pval, file=paste0("~/OneDrive/cornell/simulMI/result/gau/gene_", i, ".RData"))
    
  }
  
  return(pval)
  
}

#while(fiIdx){
#for(fiIdx in (fiIdx+1):length(file_tot)){
for (fiIdx in 1:15){
  print(fiIdx)
  print(file_tot[fiIdx])
  start_time = proc.time()
  load(paste0('/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/spot/zero_60/', file_tot[fiIdx]))
  W = matrix(0, ncol = nrow(loc), nrow =3)
  for (ct in 1:3){
    index = ((ct-1)*2500+1):(2500*ct) 
    W[ct, index] = 1
  }
  pval = tweedspot(count, W,  loc)
  end_time = proc.time()
  run_time = end_time[3] - start_time[3]
  #save(pval, gamma,run_time, file=paste0("~/OneDrive/cornell/simulMI/tweedspot/", file_tot[fiIdx] ))
  save(pval, gamma_true,run_time, file=paste0("~/OneDrive/cornell/simulMICT/bigSimul/spot/zero_60/reml/", file_tot[fiIdx] ))
}


# pval_adjusted= matrix(nrow = 100, ncol = 1)
# pval_adjusted = p.adjust(pval[1,], method = "BH")
# 
# pval_adjusted_bin <- as.matrix(pval_adjusted <0.01)
# mode(pval_adjusted_bin) <- "numeric"
# 
# metrics(gamma, pval_adjusted_bin) 

