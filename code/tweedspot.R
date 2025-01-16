library(ashr)
library(qvalue)

# transform p-values into combined p-values
file_tot = list.files('/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/reml/', '.RData')
path0 ='/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/reml/'
T_cau0 <-0
for(i in 1:length(file_tot)){
  # file-wise process
  load(paste0(path0, file_tot[i]))
  #pval[which(is.na(pval))] <- 1
  #T_cau0 = T_cau0 + tan((0.5-pval)*pi)
  #T_cau0 = T_cau0/length(file_tot)
  #P_val <- 1-pcauchy(T_cau0)
  
  G = 1000
  K = 3
  
  #pval1 = t(P_val)
  #qval1 <- matrix(qvalaue(c(pval1))$qvalue,G,K)
  pval1 = t(pval)
  pval1[which(pval1<0)]<-0
  qval1 <- matrix(qvalue(c(pval1))$qvalue,G,K)
  qval1 <- apply(pval1,2,function(x) {p.adjust(x, method = "BY")})
  
  
  #---------Cauchy combination rule---------
  
  #pval1<- matrix(NA,G,K)
  #qval1<- matrix(NA,G,K)
  #for(k in 1:K){
  #  pval1[,k] <- apply(P_val[,c(k,k+K)],1,min) # take min X/Y coordianate 
  #  qval1[,k] <- apply(Q_val[,c(k,k+K)],1,min)
  #}
  
  save(pval, pval1, qval1, run_time, gamma_true, file = paste0(path0, '/padj/', file_tot[i]))
  
}  
