
####################################   Joint Robocov model   ####################################################

library(data.table)
for(numchr in 1:22){
  df1 = data.frame(fread(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                               "Robospan_mean", "/", "100kb", "/", "100kb",
                               ".", numchr, ".annot.gz")))
  df2 = data.frame(fread(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "Robospan_mean", "/", "5kb", "/", "5kb",
                                ".", numchr, ".annot.gz")))

  df3 = data.frame(fread(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "pRobospan_mean", "/", "100kb", "/", "100kb",
                                ".", numchr, ".annot.gz")))
  df4 = data.frame(fread(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "pRobospan_mean", "/", "5kb", "/", "5kb",
                                ".", numchr, ".annot.gz")))

  df5 = data.frame(fread(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "SEG_GTEx_top10", "/", "100kb", "/", "100kb",
                                ".", numchr, ".annot.gz")))
  df6 = data.frame(fread(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "SEG_GTEx_top10", "/", "5kb", "/", "5kb",
                                ".", numchr, ".annot.gz")))

  newdf = cbind.data.frame(df1[,1:4], df1[,5], df2[,5], df3[,5], df4[,5], df5[,5], df6[,5])
  colnames(newdf) = c(colnames(df1)[1:4], "Robospan_100kb", "Robospan_5kb", "pRobospan_100kb", "pRobospan_5kb",
                      "SEG_100kb", "SEG_5kb")

  if(!dir.exists(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints_Feb2020/", "FS0"))){
    dir.create(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints_Feb2020/", "FS0"))
  }
  write.table(newdf,
              file = gzfile(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints_Feb2020/",
                                   "FS0", "/",
                                   "FS0", ".",
                                   numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  cat("We are at chrom:", numchr, "\n")


}


for(numchr in 1:22){
  df = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints_Feb2020/",
                               "FS0", "/",
                               "FS0", ".", numchr, ".annot.gz")))
  newdf = df[,c(1:4, 5, 8)]
  if(!dir.exists(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints_Feb2020/", "FS5"))){
    dir.create(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints_Feb2020/", "FS5"))
  }
  write.table(newdf,
              file = gzfile(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints_Feb2020/",
                                   "FS5", "/",
                                   "FS5", ".",
                                   numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  cat("We are at chrom:", numchr, "\n")
}

