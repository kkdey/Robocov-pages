
for(numchr in 1:22){
  df1 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "NMF_run_5_probocov", "/",
                                "NMF_run_5_clus_5", "/",
                                "100kb", "/",
                                "100kb", ".", numchr, ".annot.gz")))
  df2 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "NMF_run_5_robocov", "/",
                                "NMF_run_5_clus_3", "/",
                                "100kb", "/",
                                "100kb", ".", numchr, ".annot.gz")))
  df3 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "Robocov_hub_specific", "/",
                                "brain_specific", "/",
                                "100kb", "/",
                                "100kb", ".", numchr, ".annot.gz")))
  df4 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "Robocov_precision_hub_specific", "/",
                                "brain_specific", "/",
                                "100kb", "/",
                                "100kb", ".", numchr, ".annot.gz")))
  df5 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "SEG_GTEx_top10", "/",
                                "100kb", "/",
                                "100kb", ".", numchr, ".annot.gz")))
  newdf = cbind.data.frame(df1[,1:4], df1[,5], df2[,5], df3[,5], df4[,5], df5[,5])
  colnames(newdf) = c(colnames(df1)[1:4], "NMF_run_5_probocov_100kb",
                      "NMF_run_5_robocov_100kb", "Robocov_hub_specific_100kb",
                      "Robocov_precision_hub_specific_100kb",
                      "SEG_GTEx_top10_100kb")

  if(!dir.exists(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/Joint_Jun15/", "FS0"))){
    dir.create(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/Joint_Jun15/", "FS0"))
  }
  write.table(newdf,
              file = gzfile(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/Joint_Jun15/",
                                   "FS0", "/",
                                   "FS0", ".",
                                   numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  cat("We are at chrom:", numchr, "\n")

}

for(numchr in 1:22){
  df1 = data.frame(fread(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/Joint_Jun15/",
                                "FS0", "/",
                                "FS0", ".",
                                numchr, ".annot.gz")))
  df2 = df1[,c(1:4, 8, 9)]
  if(!dir.exists(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/Joint_Jun15/", "FS2"))){
    dir.create(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/Joint_Jun15/", "FS2"))
  }
  write.table(df2,
              file = gzfile(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/Joint_Jun15/",
                                   "FS2", "/",
                                   "FS2", ".",
                                   numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  cat("We are at chrom:", numchr, "\n")

}



df1 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                              "NMF_run_5_probocov", "/",
                              "NMF_run_5_clus_5", "/",
                              "100kb", "/",
                              "100kb", ".", numchr, ".annot.gz")))
df2 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                              "NMF_run_5_robocov", "/",
                              "NMF_run_5_clus_3", "/",
                              "100kb", "/",
                              "100kb", ".", numchr, ".annot.gz")))
df3 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                              "Robocov_hub_specific", "/",
                              "brain_specific", "/",
                              "100kb", "/",
                              "100kb", ".", numchr, ".annot.gz")))
df4 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                              "Robocov_precision_hub_specific", "/",
                              "brain_specific", "/",
                              "100kb", "/",
                              "100kb", ".", numchr, ".annot.gz")))
df5 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                              "Cor_hub_specific", "/",
                              "brain_specific", "/",
                              "100kb", "/",
                              "100kb", ".", numchr, ".annot.gz")))
df6 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                              "NMF_run_5_standard", "/",
                              "NMF_run_5_clus_5", "/",
                              "100kb", "/",
                              "100kb", ".", numchr, ".annot.gz")))

newdf = cbind(df1[,5], df2[,5], df3[,5], df4[,5], df5[,5], df6[,5])








