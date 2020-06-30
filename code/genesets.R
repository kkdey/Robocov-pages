

df = read.table("/Users/kushaldey/Documents/Robocov-pages/output/gene_anno_unique_datefix.txt", header=T)

filename = "NMF_run_5_probocov"
ll = list.files(paste0("/Users/kushaldey/Documents/Robocov-pages/output/", filename), pattern = ".txt")
for(num in 1:length(ll)){
  test1 = as.character(read.delim(paste0("/Users/kushaldey/Documents/Robocov-pages/output/",
                                 filename, "/",
                                 ll[num]),
                                  header = F)[,1])
  tmp = cbind.data.frame(df$symbol[match(intersect(df$id, test1), df$id)], 1)
  if(!dir.exists(paste0("/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/", filename))){
    dir.create(paste0("/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/", filename))
  }
  write.table(tmp, file = paste0("/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/",
                          filename, "/",
                          ll[num]),
              row.names = F, col.names = F, sep = "\t", quote=F)

  cat("We read filename:", numl)
}


filename = "Cor_hub_specific"
ll = list.files(paste0("/Users/kushaldey/Documents/Robocov-pages/output/", filename), pattern = ".txt")
ll= ll[-c(8,9)]
for(num in 1:length(ll)){
  test1 = as.character(read.delim(paste0("/Users/kushaldey/Documents/Robocov-pages/output/",
                                         filename, "/",
                                         ll[num]),
                                  header = F)[,1])
  tmp = cbind.data.frame(df$symbol[match(intersect(df$id, test1), df$id)], 1)
  if(!dir.exists(paste0("/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/", filename))){
    dir.create(paste0("/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/", filename))
  }
  write.table(tmp, file = paste0("/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/",
                                 filename, "/",
                                 ll[num]),
              row.names = F, col.names = F, sep = "\t", quote=F)

  cat("We read filename:", numl)
}





test1 = as.character(read.delim("/Users/kushaldey/Documents/Robocov-pages/output/NMF_run_5_robocov/NMF_run_5_clus_5.txt",
                   header = F)[,1])
test2 = as.character(read.delim("/Users/kushaldey/Documents/Robocov-pages/output/NMF_run_5_probocov/NMF_run_5_clus_5.txt",
                                header = F)[,1])
test3 = as.character(read.delim("/Users/kushaldey/Documents/Robocov-pages/output/NMF_run_5_standard/NMF_run_5_clus_4.txt",
                                header = F)[,1])

tmp = cbind.data.frame(df$symbol[match(intersect(df$id, test1), df$id)], 1)
write.table(tmp, file = "/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/NMF_run_5_robocov_clust_5.txt",
            row.names = F, col.names = F, sep = "\t", quote=F)


tmp = cbind.data.frame(df$symbol[match(intersect(df$id, test2), df$id)], 1)
write.table(tmp, file = "/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/NMF_run_5_probocov_clust_5.txt",
            row.names = F, col.names = F, sep = "\t", quote=F)

tmp = cbind.data.frame(df$symbol[match(intersect(df$id, test3), df$id)], 1)
write.table(tmp, file = "/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/NMF_run_5_standard_clust_4.txt",
            row.names = F, col.names = F, sep = "\t", quote=F)





test1 = as.character(read.delim("/Users/kushaldey/Documents/Robocov-pages/output/Cor_hub_specific/brain_specific.txt",
                                header = F)[,1])
test2 = as.character(read.delim("/Users/kushaldey/Documents/Robocov-pages/output/Robocov_hub_specific/brain_specific.txt",
                                header = F)[,1])
test3 = as.character(read.delim("/Users/kushaldey/Documents/Robocov-pages/output/Robocov_precision_hub_specific/brain_specific.txt",
                                header = F)[,1])

tmp = cbind.data.frame(df$symbol[match(intersect(df$id, test1), df$id)], 1)
write.table(tmp, file = "/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/Hub_brain_cor.txt",
            row.names = F, col.names = F, sep = "\t", quote=F)


tmp = cbind.data.frame(df$symbol[match(intersect(df$id, test2), df$id)], 1)
write.table(tmp, file = "/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/Hub_brain_robocov.txt",
            row.names = F, col.names = F, sep = "\t", quote=F)

tmp = cbind.data.frame(df$symbol[match(intersect(df$id, test3), df$id)], 1)
write.table(tmp, file = "/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/Hub_brain_probocov.txt",
            row.names = F, col.names = F, sep = "\t", quote=F)
