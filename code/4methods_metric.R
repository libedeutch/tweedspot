# TDR, FDR v.s. Level alpha set-up 

library(ggplot2)
library(pROC)
metrics <- function(gamma_true, pval_adjusted_bin){
  
  tp = sum( gamma_true==1 & pval_adjusted_bin==1 )
  tn = sum( gamma_true==0 & pval_adjusted_bin==0 )
  fp = sum( gamma_true==0 & pval_adjusted_bin==1 )
  fn = sum( gamma_true==1 & pval_adjusted_bin==0 )
  
  power = 1-fn/(fn+tp+0.00000001)
  fdr = fp/(fp+tp+0.00000001)
  return(list(tp=tp, tn = tn, fp = fp, fn = fn, power = power, fdr =fdr))
}
library(stringr)


#qval1 <- apply(pval,1,function(x) {p.adjust(x, method = "BY")})

cut_offs = c(0.01,0.05,0.10,0.20,0.25)
perf_metrics = c()
for (cut_off in cut_offs){
  
  path = "/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/sparkx/"
  file_tot = list.files(path, ".RData")
  perf_metrics_sparkx = matrix(NA, nrow = length(file_tot), ncol = 13)
  colnames(perf_metrics_sparkx) = c("TP", "TN", "FP", "FN", "power", "FDR", "pattern", "zero", "replicate","method", 'run_time', 'AUC','alpha')
  for(i in 1:length(file_tot)){
    load(paste0(path,file_tot[i]))
    #load(paste0(path_old, file_tot[i]))
    #pval_adjusted =apply(pval,1,function(x) {p.adjust(x, method = "BY")})   
    #pval_adjusted_bin <- as.matrix(pval_adjusted <0.01)
    #pval_adjusted_bin <- p.adjust(CombinePValues(t(pval)), method = 'BY')<0.01
    #pval_adjusted_bin <- p.adjust(apply(pval,2,min))<0.01
    #mode(pval_adjusted_bin) <- "numeric"
    perf_metrics_sparkx[i, 1:6] =  unlist(metrics(apply(gamma_true,2,max), t(result$combinedPval<cut_off) ))
    #perf_metrics_tweedspot[i, 1:6] =  unlist(metrics(apply(gamma_true,2,max), pval_adjusted_bin ))
    perf_metrics_sparkx[i,7] = substr(file_tot[i], 1, str_locate(file_tot[i], "pattern")[1]-2 )
    perf_metrics_sparkx[i,8] = substr(file_tot[i], str_locate(file_tot[i], "zero_")[2]+1, str_locate(file_tot[i], "_replicate")[1]-1 )
    perf_metrics_sparkx[i,9] = substr(file_tot[i], str_locate(file_tot[i], "replicate_")[2]+1, str_locate(file_tot[i], ".RData")[1]-1 )
    perf_metrics_sparkx[i,10]  ="SPARKX"
    perf_metrics_sparkx[i,11]  =ifelse(is.null(run_time), 0,run_time)
    roc_curve <- roc(response = as.vector(apply(gamma_true,2,max)), predictor = as.vector(t(ifelse(result$combinedPval<cut_off,1,0) )))
    # Calculate the AUC
    perf_metrics_sparkx[i,12]  = auc(roc_curve)
    perf_metrics_sparkx[i,13]  = cut_off
  }
  
  path = "/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/sparkx/sparkx_ct/"
  file_tot = list.files(path, ".RData")
  perf_metrics_sparkx_ct = matrix(NA, nrow = length(file_tot), ncol = 13)
  colnames(perf_metrics_sparkx_ct) = c("TP", "TN", "FP", "FN", "power", "FDR", "pattern", "zero", "replicate","method", 'run_time', 'AUC','alpha')
  for(i in 1:length(file_tot)){
    load(paste0(path,file_tot[i]))
    #load(paste0(path_old, file_tot[i]))
    #pval_adjusted =apply(pval,1,function(x) {p.adjust(x, method = "BY")})
    #pval_adjusted_bin <- as.matrix(pval_adjusted <0.01)
    #pval_adjusted_bin <- p.adjust(CombinePValues(t(pval)), method = 'BY')<0.01
    #pval_adjusted_bin <- p.adjust(apply(pval,2,min))<0.01
    #mode(pval_adjusted_bin) <- "numeric"
    perf_metrics_sparkx_ct[i, 1:6] =  unlist(metrics(apply(gamma_true,2,max), t(result$combinedPval<cut_off) ))
    #perf_metrics_tweedspot[i, 1:6] =  unlist(metrics(apply(gamma_true,2,max), pval_adjusted_bin ))
    perf_metrics_sparkx_ct[i,7] = substr(file_tot[i], 1, str_locate(file_tot[i], "pattern")[1]-2 )
    perf_metrics_sparkx_ct[i,8] = substr(file_tot[i], str_locate(file_tot[i], "zero_")[2]+1, str_locate(file_tot[i], "_replicate")[1]-1 )
    perf_metrics_sparkx_ct[i,9] = substr(file_tot[i], str_locate(file_tot[i], "replicate_")[2]+1, str_locate(file_tot[i], ".RData")[1]-1 )
    perf_metrics_sparkx_ct[i,10]  ="SPARKXCT"
    perf_metrics_sparkx_ct[i,11]  =ifelse(is.null(run_time), 0,run_time)
    roc_curve <- roc(response = as.vector(apply(gamma_true,2,max)), predictor = as.vector(t(ifelse(result$combinedPval<cut_off,1,0) )))
    # Calculate the AUC
    perf_metrics_sparkx_ct[i,12]  = auc(roc_curve)
    perf_metrics_sparkx_ct[i,13]  = cut_off
  }
  
  path = "/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/nnsvg/"
  file_tot = list.files(path, ".RData")
  perf_metrics_nnsvg = matrix(NA, nrow = length(file_tot), ncol = 13)
  colnames(perf_metrics_nnsvg) = c("TP", "TN", "FP", "FN", "power", "FDR", "pattern", "zero", "replicate","method", 'run_time', 'AUC','alpha')
  for(i in 1:length(file_tot)){
    load(paste0(path,file_tot[i]))
    #load(paste0(path_old, file_tot[i]))
    #pval_adjusted =apply(pval,1,function(x) {p.adjust(x, method = "BY")})
    #pval_adjusted_bin <- as.matrix(pval_adjusted <0.01)
    #pval_adjusted_bin <- p.adjust(CombinePValues(t(pval)), method = 'BY')<0.01
    #pval_adjusted_bin <- p.adjust(apply(pval,2,min))<0.01
    #mode(pval_adjusted_bin) <- "numeric"
    perf_metrics_nnsvg[i, 1:6] =  unlist(metrics(apply(gamma_true,2,max), t(cc[,'padj']<cut_off) ))
    #perf_metrics_tweedspot[i, 1:6] =  unlist(metrics(apply(gamma_true,2,max), pval_adjusted_bin ))
    perf_metrics_nnsvg[i,7] = substr(file_tot[i], 1, str_locate(file_tot[i], "pattern")[1]-2 )
    perf_metrics_nnsvg[i,8] = substr(file_tot[i], str_locate(file_tot[i], "zero_")[2]+1, str_locate(file_tot[i], "_replicate")[1]-1 )
    perf_metrics_nnsvg[i,9] = substr(file_tot[i], str_locate(file_tot[i], "replicate_")[2]+1, str_locate(file_tot[i], ".RData")[1]-1 )
    perf_metrics_nnsvg[i,10]  ="nnSVG"
    perf_metrics_nnsvg[i,11]  =ifelse(is.null(run_time), 0,run_time)
    roc_curve <- roc(response = as.vector(apply(gamma_true,2,max)), predictor = as.vector(t(ifelse(cc[,'padj']<cut_off,1,0) )))
    # Calculate the AUC
    perf_metrics_nnsvg[i,12]  = auc(roc_curve)
    perf_metrics_nnsvg[i,13]  = cut_off
  }
  
  
  path = "/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/reml/padj/"
  file_tot = list.files(path, ".RData")
  perf_metrics_tweedspot = matrix(NA, nrow = length(file_tot), ncol = 13)
  colnames(perf_metrics_tweedspot) = c("TP", "TN", "FP", "FN", "power", "FDR", "pattern", "zero", "replicate","method", 'run_time', 'AUC','alpha')
  for(i in 1:length(file_tot)){
    load(paste0(path,file_tot[i]))
    #load(paste0(path_old, file_tot[i]))
    #pval_adjusted =apply(pval,1,function(x) {p.adjust(x, method = "BY")})   
    #pval_adjusted_bin <- as.matrix(pval_adjusted <0.01)
    #pval_adjusted_bin <- p.adjust(CombinePValues(t(pval)), method = 'BY')<0.01
    #pval_adjusted_bin <- p.adjust(apply(pval,2,min))<0.01
    #mode(pval_adjusted_bin) <- "numeric"
    perf_metrics_tweedspot[i, 1:6] =  unlist(metrics(gamma_true, t(qval1<cut_off) ))
    #perf_metrics_tweedspot[i, 1:6] = unlist(metrics(apply(gamma_true,2,max), apply(t(qval1<cut_off),2,max)))
    #perf_metrics_tweedspot[i, 1:6] =  unlist(metrics(apply(gamma_true,2,max), pval_adjusted_bin ))
    perf_metrics_tweedspot[i,7] = substr(file_tot[i], 1, str_locate(file_tot[i], "pattern")[1]-2 )
    perf_metrics_tweedspot[i,8] = substr(file_tot[i], str_locate(file_tot[i], "zero_")[2]+1, str_locate(file_tot[i], "_replicate")[1]-1 )
    perf_metrics_tweedspot[i,9] = substr(file_tot[i], str_locate(file_tot[i], "replicate_")[2]+1, str_locate(file_tot[i], ".RData")[1]-1 )
    perf_metrics_tweedspot[i,10]  ="TweedSpot"
    perf_metrics_tweedspot[i,11]  =ifelse(is.null(run_time), 0,run_time)
    roc_curve <- roc(response = as.vector(gamma_true), predictor = as.vector(t(ifelse(qval1<cut_off,1,0) )))
    #roc_curve <- roc(response = apply(gamma_true,2,max), predictor = apply(t(ifelse(qval1<cut_off,1,0) ),2,max))
    
    # Calculate the AUC
    perf_metrics_tweedspot[i,12]  = auc(roc_curve)
    perf_metrics_tweedspot[i,13]  = cut_off
  }
  
  path = "/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/ctsv/zero_30/padj2/"
  file_tot = list.files(path, ".RData")
  perf_metrics_ctsv = matrix(NA, nrow = length(file_tot), ncol = 13)
  colnames(perf_metrics_ctsv) = c("TP", "TN", "FP", "FN", "power", "FDR", "pattern", "zero", "replicate","method", 'run_time', 'AUC','alpha')
  for(i in 1:length(file_tot)){
    load(paste0(path,file_tot[i]))
    #load(paste0(path_old, file_tot[i]))
    #pval_adjusted =apply(pval,1,function(x) {p.adjust(x, method = "BY")})   
    #pval_adjusted_bin <- as.matrix(pval_adjusted <0.01)
    #pval_adjusted_bin <- p.adjust(CombinePValues(t(pval)), method = 'BY')<0.01
    #pval_adjusted_bin <- p.adjust(apply(pval,2,min))<0.01
    #mode(pval_adjusted_bin) <- "numeric"
    perf_metrics_ctsv[i, 1:6] =  unlist(metrics(gamma_true, t(qval1<cut_off) ))
    #perf_metrics_tweedspot[i, 1:6] =  unlist(metrics(apply(gamma_true,2,max), pval_adjusted_bin ))
    perf_metrics_ctsv[i,7] = substr(file_tot[i], 1, str_locate(file_tot[i], "pattern")[1]-2 )
    perf_metrics_ctsv[i,8] = substr(file_tot[i], str_locate(file_tot[i], "zero_")[2]+1, str_locate(file_tot[i], "_replicate")[1]-1 )
    perf_metrics_ctsv[i,9] = substr(file_tot[i], str_locate(file_tot[i], "replicate_")[2]+1, str_locate(file_tot[i], ".RData")[1]-1 )
    perf_metrics_ctsv[i,10]  ="CTSV"
    perf_metrics_ctsv[i,11]  =ifelse(is.null(run_time), 0,run_time)
    roc_curve <- roc(response = as.vector(gamma_true), predictor = as.vector(t(ifelse(qval1<cut_off,1,0) )))
    # Calculate the AUC
    perf_metrics_ctsv[i,12]  = auc(roc_curve)
    perf_metrics_ctsv[i,13]  = cut_off
  }
  
  
  #perf_metrics = rbind(perf_metrics, perf_metrics_sparkx, perf_metrics_sparkx_ct, perf_metrics_nnsvg, perf_metrics_tweedspot, perf_metrics_ctsv)
  perf_metrics = rbind(perf_metrics, perf_metrics_sparkx, perf_metrics_tweedspot, perf_metrics_ctsv)
}



#perf_metrics = rbind(perf_metrics_spkx, perf_metrics_nnsvg, perf_metrics_tweedspot)
#perf_metrics = rbind(perf_metrics_ctsv, perf_metrics_tweedspot)

library(ggplot2)

plotdf = data.frame(perf_metrics)


#pattern = 'linear' # [which(plotdf$pattern==pattern ),]
power_plot_BY = ggplot(plotdf,aes(x = as.factor(method),y = as.numeric(power)))+
  geom_boxplot(aes(fill = as.factor(method)))+
  facet_grid(as.factor(pattern)~as.factor(alpha),)+
  #xlab(paste0(pattern, " Pattern"))+
  ylab("Power")+
  theme(
    strip.text = element_text(size = 14),
    strip.background =element_rect(fill = 'lightgreen'),
    #panel.background = element_blank(),
    legend.position = "top")+ggtitle('Power')+labs(fill = "Method") 

fdr_plot_BY = ggplot(plotdf,aes(x = as.factor(method),y = as.numeric(FDR)))+
  geom_boxplot(aes(fill = as.factor(method)))+
  geom_hline(data = plotdf, aes(yintercept = as.numeric(alpha)), color = "red", linetype = "dashed", size = 1) +
  #geom_hline(yintercept = as.factor(alpha), color = "red", linetype = "dashed", size = 1) + 
  facet_grid(as.factor(pattern)~~as.factor(alpha),)+
  #xlab(paste0(pattern, " Pattern"))+
  ylab("fdr")+
  theme(
    strip.text = element_text(size = 14),
    strip.background =element_rect(fill = 'lightgreen'),
    #panel.background = element_blank(),
    legend.position = "top") + ggtitle("FDR") +labs(fill = "Method") 


auc_plot_BY = ggplot(plotdf[which(plotdf$alpha==0.05),],aes(x = as.factor(method),y = as.numeric(AUC)))+
  geom_violin(aes(fill = as.factor(method)))+
  #geom_hline(data = plotdf, aes(yintercept = as.numeric(alpha)), color = "red", linetype = "dashed", size = 1) +
  #geom_hline(yintercept = as.factor(alpha), color = "red", linetype = "dashed", size = 1) + 
  #facet_grid(as.factor(alpha)~as.factor(pattern),)+
  facet_grid(as.factor(pattern)~as.factor(alpha),)+
  #xlab(paste0(pattern, " Pattern"))+
  ylab("AUC")+
  theme( legend.position = 'top',
         strip.text = element_text(size = 14),
         strip.background =element_rect(fill = 'lightgreen')
         #panel.background = element_blank(),
  ) + ggtitle("AUC")+labs(fill = "Method") 
cowplot::plot_grid(power_plot_BY,  fdr_plot_BY)#auc_plot_BY,


fdr_plot_alpha = ggplot(plotdf,aes(x = as.numeric(alpha) ,y = as.numeric(FDR)))+
  geom_line(aes(x = as.numeric(alpha) ,y = as.numeric(FDR)))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +
  facet_grid(as.factor(method)~as.factor(pattern),)+
  xlab(" ")+
  ylab("FDR")+
  theme(strip.text = element_text(size = 14),
        strip.background =element_rect(fill = 'lightgreen')
        #panel.background = element_blank(),
  )+ggtitle('FDR')#+labs(fill = "Method") 

tdr_plot_alpha = ggplot(plotdf,aes(x = as.numeric(alpha) ,y = 1-as.numeric(FDR)))+
  geom_line(aes(x = as.numeric(alpha) ,y = 1-as.numeric(FDR)))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +
  facet_grid(as.factor(method)~as.factor(pattern),)+
  xlab(" ")+
  ylab("TDR")+
  theme(strip.text = element_text(size = 14),
        strip.background =element_rect(fill = 'lightgreen'),
        #panel.background = element_blank(),
  )+ggtitle('TDR')

library(dplyr)
plotdf2 = plotdf %>% group_by(method, alpha)%>%
  summarize(mean_fdr = mean(as.numeric(FDR)), mean_power = mean(as.numeric(power)), mean_auc = mean(as.numeric(AUC)))

ggplot(plotdf2,aes(x = as.numeric(mean_fdr) ,y = as.numeric(mean_auc), color = as.factor(method)))+
  geom_line(size = 1)+
  #geom_point()+
  #geom_boxplot()+
  #geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +
  #facet_grid(as.factor(pattern)~.,)+
  xlab(" ")+
  ylab("power")+
  theme(strip.text = element_text(size = 14),
        strip.background =element_rect(fill = 'lightgreen')
        #panel.background = element_blank(),
  )+ggtitle('FDR')

