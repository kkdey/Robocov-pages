

library(data.table)

folder = paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/SEG_GTEx_top10")
ll = list.files(folder)


total_vec = c()
for(mm in 1:length(ll)){
  vec = c()
  for(numchr in 1:22){
    df1 = data.frame(fread(paste0("zcat ", folder, "/",
                                  ll[mm], "/",
                                  ll[mm], ".", numchr, ".annot.gz")))
    vec = c(vec, df1[,5])
  }
  total_vec = c(total_vec, sum(vec)/length(vec))
}

names(total_vec) = ll
total_vec*100
