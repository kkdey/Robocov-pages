
for(numchr in 1:22){
  df1 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "Robospan_mean", "/",
                                "100kb", "/",
                                "100kb", ".", numchr, ".annot.gz")))
  df2 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "Robospan_mean", "/",
                                "5kb", "/",
                                "5kb", ".", numchr, ".annot.gz")))
  df3 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "Robospan_mean", "/",
                                "CODING", "/",
                                "CODING", ".", numchr, ".annot.gz")))
  df4 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "Robospan_mean", "/",
                                "TSS", "/",
                                "TSS", ".", numchr, ".annot.gz")))
  df5 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "pRobospan_mean", "/",
                                "100kb", "/",
                                "100kb", ".", numchr, ".annot.gz")))
  df6 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "pRobospan_mean", "/",
                                "5kb", "/",
                                "5kb", ".", numchr, ".annot.gz")))
  df7 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "pRobospan_mean", "/",
                                "CODING", "/",
                                "CODING", ".", numchr, ".annot.gz")))
  df8 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "pRobospan_mean", "/",
                                "Roadmap_Enhancer", "/",
                                "Roadmap_Enhancer", ".", numchr, ".annot.gz")))
  df9 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "pRobospan_mean", "/",
                                "TSS", "/",
                                "TSS", ".", numchr, ".annot.gz")))

  newdf = cbind.data.frame(df1[,1:4], df1[,5], df2[,5], df3[,5], df4[,5], df5[,5], df6[,5], df7[,5], df8[,5], df9[,5])
  colnames(newdf) = c(colnames(df1)[1:4], "Robo_100", "Robo_5", "Robo_Coding", "Robo_TSS",
                      "pRobo_100", "pRobo_5", "pRobo_Coding", "pRobo_Road", "pRobo_TSS")

  if(!dir.exists(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/Joint_Nov29/", "FS0"))){
    dir.create(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/Joint_Nov29/", "FS0"))
  }
  write.table(newdf,
              file = gzfile(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/Joint_Nov29/",
                                   "FS0", "/",
                                   "FS0", ".",
                                   numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  cat("We are at chrom:", numchr, "\n")

}


for(numchr in 1:22){
  df = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/Joint_Nov29/",
                                "FS1", "/",
                                "FS1", ".", numchr, ".annot.gz")))
  newdf = df[, -c(6,7)]
  if(!dir.exists(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/Joint_Nov29/", "FS3"))){
    dir.create(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/Joint_Nov29/", "FS3"))
  }
  write.table(newdf,
              file = gzfile(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/Joint_Nov29/",
                                   "FS3", "/",
                                   "FS3", ".",
                                   numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  cat("We are at chrom:", numchr, "\n")

}
