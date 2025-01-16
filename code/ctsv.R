# CTSV 

library(DIRECT)
library(pscl)
#library(doSNOW)
#library(qvalue)

###############################################################
########################### CTSV ##############################
###############################################################
#--------
#load("/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/linear/linear_pattern_zero_0_replicate_1.RData")
# normalize coordinates

file_tot = list.files('/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/spot/zero_60/', ".RData")
path_d = '~/OneDrive/cornell/simulMICT/bigSimul/spot/zero_60/'
path0 = '~/OneDrive/cornell/simulMICT/bigSimul/spot/zero_60/ctsv/'
ctsv <- function(Y,W,loc  ){
  S <- t(loc) - colMeans(loc)
  S <- t(S / apply(S, 1, sd))
  
  quan <- c(0.2,0.4,0.6,0.8,1.0)
  psi1 <- quantile(abs(S[,1]), quan)
  psi2 <- quantile(abs(S[,2]), quan)
  #spot number
  n <- nrow(loc)

  
  #gene number
  G <- nrow(Y)
  
  # cell type number
  K <- 3
  
  
  for(fit_pat in c("gau2","cos2","gau1","cos1","linear")){
    if(fit_pat == "gau1"){
      h1 <- exp(-S[,1]^2 / 2 / psi1[2]^2)
      h2 <- exp(-S[,2]^2 / 2 / psi2[2]^2)
    }else if(fit_pat == "gau2"){
      h1 <- exp(-S[,1]^2 / 2 / psi1[3]^2)
      h2 <- exp(-S[,2]^2 / 2 / psi2[3]^2)
    }else if(fit_pat == "cos1"){
      h1 <- cos(2*pi*S[,1] / psi1[2])
      h2 <- cos(2*pi*S[,2] / psi2[2])
    } else if(fit_pat == "cos2"){
      h1 <- cos(2*pi*S[,1] / psi1[3])
      h2 <- cos(2*pi*S[,2] / psi2[3])
    }else{
      h1 <- S[,1]
      h2 <- S[,2]
    }
    
    Tmp <- cbind(t(W) * h1, t(W) * h2, t(W))
    colnames(Tmp) <- 1:ncol(Tmp)
    start_time = proc.time()
    pval = matrix(nrow = 2*K, ncol = G)
    for(g in 1:nrow(Y)){
      y <- Y[g,]
      ell <- colSums(Y) / median(colSums(Y))
      fm_zinb0 <- zeroinfl(y ~ -1+offset(log(ell))+Tmp|1,
                           dist = "negbin",link = "probit",
                           control = zeroinfl.control(method = "CG"
                           ))
      p_val <- coef(summary(fm_zinb0))$count[,4]
      nind <- 2*(length(p_val) - 1)/3
      p_val <- p_val[1:nind]
      pval[1:(2*K),g] = p_val
  }
    end_time = proc.time()
    run_time = end_time[3] - start_time[3]
    save(pval,run_time, gamma_true, file = paste0(path0, fit_pat, '_', file_tot[fiIdx]) )
  
  
  }
  
} 

for (fiIdx in 1:15){
  print(fiIdx)
  print(file_tot[fiIdx])
  #start_time = proc.time()
  load(paste0(path_d, file_tot[fiIdx]))
  W = matrix(0, ncol = nrow(loc), nrow =3)
  for (ct in 1:3){
    index = ((ct-1)*2500+1):(2500*ct)
    W[ct, index] = 1
  }
  ctsv(t(count), W,  loc)
  #end_time = proc.time()
  #run_time = end_time[3] - start_time[3]
  #save(pval, gamma,run_time, file=paste0("~/OneDrive/cornell/simulMI/tweedspot/", file_tot[fiIdx] ))
  #save(pval, gamma_true,run_time, file=paste0("~/OneDrive/cornell/simulMICT/ctsv/", file_tot[fiIdx] ))
}











