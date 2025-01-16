install.packages("BiocManager")
BiocManager::install("nnSVG")
install.packages("BiocManager")
BiocManager::install("SpatialExperiment")
BiocManager::install("STexampleData")
install.packages("BRISC")

library(nnSVG)
library(STexampleData)
library(scran)
library(ggplot2)



path = '/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/spot/zero_60/'
save_path = '/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul/spot/zero_60/nnsvg/'
file_tot = list.files(path, '.RData')
idx = c(1:3, 4:length(file_tot))
for(i in idx){
  # nnsvg
  print(i)
  start_time = proc.time()  
  load(paste0(path, file_tot[i] ))
  logNcount = t(count)
  logNcount = apply(logNcount,2,function(x){x/sum(x)})
  logNcount = log(logNcount+1)
  cc = nnSVG(input = logNcount, spatial_coords = loc )
  end_time = proc.time()
  run_time = end_time[3] - start_time[3]
  save(cc,run_time,gamma_true , file = paste0(save_path,file_tot[i]) )
}
