file_tot = list.files('/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/spot/zero_60/reml', '.RData')
path = '/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/spot/zero_60/ctsv/'

for (i in 1:length(file_tot)){
  runtime = 0
  load(paste0(path, 'gau2_', file_tot[i]))
  pval[which(is.na(pval))] <- 1
  p_gau2 = t(pval)
  runtime = runtime+run_time
  rm(pval)
  rm(run_time)
  
  load(paste0(path, 'gau1_', file_tot[i]))
  pval[which(is.na(pval))] <- 1
  p_gau1 = t(pval)
  runtime = runtime+run_time
  rm(pval)
  rm(run_time)
  
  load(paste0(path, 'cos2_', file_tot[i]))
  pval[which(is.na(pval))] <- 1
  p_cos2 = t(pval)
  runtime = runtime+run_time
  rm(pval)
  rm(run_time)
  
  load(paste0(path, 'cos1_', file_tot[i]))
  pval[which(is.na(pval))] <- 1
  p_cos1 = t(pval)
  runtime = runtime+run_time
  rm(pval)
  rm(run_time)
  
  load(paste0(path, 'linear_', file_tot[i]))
  pval[which(is.na(pval))] <- 1
  p_linear = t(pval)
  runtime = runtime+run_time
  
  G = ncol(pval)
  K = nrow(pval)/2
  rm(pval)
  rm(run_time)
  
  
  
  #---------Cauchy combination rule---------
  
  T_cau0 <- (tan((0.5-p_gau2) * pi) + tan((0.5-p_gau1) * pi) + 
               tan((0.5-p_cos2) * pi)+
               tan((0.5-p_cos1) * pi) + 
               tan((0.5-p_linear) * pi) ) / 5
  P_val_ZINB <- 1-pcauchy(T_cau0)
  
  if(FALSE){
  Q_ZINB1 <- matrix(qvalue(c(P_val_ZINB))$qvalue,G,2*K)
  
  P_val_ZINB1<- matrix(NA,G,K)
  Q_val_ZINB1<- matrix(NA,G,K)
  for(k in 1:K){
    P_val_ZINB1[,k] <- apply(P_val_ZINB[,c(k,k+K)],1,min)
    Q_val_ZINB1[,k] <- apply(Q_ZINB1[,c(k,k+K)],1,min)
  }
  }
  run_time = runtime
  pval = P_val_ZINB

  save(pval, gamma_true, run_time, file = paste0(path,'padj_base/', file_tot[i]))
}