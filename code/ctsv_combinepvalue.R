library(qvalue)
library(SPARK)

# transform p-values into combined p-values
file_tot = list.files('/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/spot/zero_60/ctsv/padj_base/', '.RData')
path0 ='/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/spot/zero_60/ctsv/padj_base/'
path1 ='/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/spot/zero_60/ctsv/padj/'

for(i in 1:length(file_tot)){
  # file-wise process
  load(paste0(path0, file_tot[i]))

  G = 1000
  K = 3
  
  #P_val = t(P_val)
  #Q_val <- matrix(qvalue(c(P_val))$qvalue,G,2*K)
  
  P_val = pval 
  #Q_val <- matrix(qvalue(c(P_val))$qvalue,G,2*K)
  Q_val <- apply(P_val,2,function(x) {p.adjust(x, method = "BY")})  
  #Q_val <- p.adjust(P_val, method = 'BY')
  #---------Cauchy combination rule---------
  
  pval1<- matrix(NA,G,K)
  qval1<- matrix(NA,G,K)
  for(k in 1:K){
    pval1[,k] <- apply(P_val[,c(k,k+K)],1,min) # take min X/Y coordianate 
    qval1[,k] <- apply(Q_val[,c(k,k+K)],1,min)
  }

save(pval, pval1, qval1, run_time, gamma_true, file = paste0(path1,  file_tot[i]))
  
}  



