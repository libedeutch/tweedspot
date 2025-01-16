library(SPARK)
path = '/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/spot/zero_60/'
save_path = '/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/spot/zero_60/sparkx/sparkx_ct/'
file_tot = list.files(path, '.RData')
for(i in 1:length(file_tot)){
  # nnsvg
  start_time = proc.time()
  load(paste0(path, file_tot[i] ))
  S <- t(loc) - colMeans(loc)
  S <- t(S / apply(S, 1, sd))
  for (ct in 1:3){
    index = ((ct-1)*2500+1):(2500*ct)
    #W[ct, index] = 1
    count_sub = t(count[index,])
    rownames(count_sub) <- paste0('spot_', 1:nrow(count_sub))
    S_sub = S[index,]
    cc = sparkx(count_sub,S_sub)
    end_time = proc.time()
    run_time = end_time - start_time
    result = cc$res_mtest
    run_time = run_time[3]
    save(result,gamma_true, run_time, file = paste0(save_path,'ct_',ct,'_',file_tot[i]) )
  }

}

library(stringr)
Path = '/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/spot/zero_60/sparkx/'
Path2 = '/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/spot/zero_60/sparkx/sparkx_ct/'
file_tot = list.files(Path, '.RData')

for(i in 1:length(file_tot)){
  pval = c()
  run_time_s = 0
  for(j in 1:3){
    load(paste0(Path2, 'ct_', j, '_', file_tot[i]))
    
    spots = as.vector(as.numeric(str_remove_all(rownames(result), 'spot_')))
    missing_spot = c(1:1000)[!(c(1:1000) %in% spots)]
    
    result0 = matrix(ncol = 2, nrow = 1000)
    result0[missing_spot,] = c(1,1)
    result0[spots, ] = as.matrix(result) #[1:nrow(result),]
    
    colnames(result0) = colnames(result)
    result = data.frame(result0)
    
    pval = rbind(pval,result$combinedPval)
    run_time_s = run_time_s + run_time
  }
  run_time = run_time_s
  save(pval, run_time, gamma_true,file = paste0(Path2, 'padj/', file_tot[i]))
}

