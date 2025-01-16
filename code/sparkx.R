library(SPARK)
path = '/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/spot/zero_60/'
save_path = '/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/spot/zero_60/sparkx/'
file_tot = list.files(path, '.RData')
for(i in 1:length(file_tot)){
  # nnsvg
  start_time = proc.time()
  load(paste0(path, file_tot[i] ))
  cc = sparkx(t(count),loc)
  end_time = proc.time()
  run_time = end_time - start_time
  result = cc$res_mtest
  run_time = run_time[3]
  save(result,gamma_true, run_time, file = paste0(save_path,file_tot[i]) )
}

sum(is.na(count))
sum(is.nan(count))
sum(is.infinite(count))


# TP FP TN FN 
metrics <- function(gamma_true, pval_adjusted_bin){
  
  tp = sum( gamma_true==1 & pval_adjusted_bin==1 )
  tn = sum( gamma_true==0 & pval_adjusted_bin==0 )
  fp = sum( gamma_true==0 & pval_adjusted_bin==1 )
  fn = sum( gamma_true==1 & pval_adjusted_bin==0 )
  
  power = 1-fn/(fn+tp)
  fdr = fp/(fp+tp)
  return(list(tp=tp, tn = tn, fp = fp, fn = fn, power = power, fdr =fdr))
}

metrics(apply(gamma_true, 2 ,max), (result$combinedPval<0.01))

metrics(gamma_true[1,], (result$combinedPval<0.01))
