library(rmeta)
foldername="/Users/kushaldey/Documents/COVID19/data/MAGMA/magma_genesets/robocov_nmf/"
ll = list.files(foldername, pattern = ".gsa.out")
newdf_beta = c()
newdf_beta_sd = c()
library(qvalue)

outlist = vector(mode = "list", length(ll))
for(mm in 1:length(ll)){
  temp = read.table(paste0(foldername, "/", ll[mm]), header = T)
  temp = temp[order(temp$P, decreasing = F), ]
  qvalues_temp = qvalue(temp$P, pi0=1)$qvalues
  #idx = which(qvalues_temp < 0.05)
  idx = which(qvalues_temp < 0.05)
  outlist[[mm]] = as.character(temp$FULL_NAME)[idx]
}

names(outlist) = ll
len = unlist(lapply(outlist, function(x) length(x)))
outlist2 = outlist[which(len >= 1)]
outlist3 = outlist2[-grep("UKB2_145K", names(outlist2))]

tt = as.character(unlist(outlist3))
length(grep("probocov", tt))
length(grep("robocov", tt))
length(grep("standard", tt))
length(grep("cor", tt))
