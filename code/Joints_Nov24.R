
library(data.table)

for(numchr in 1:22){
  df1 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "Robospan_mean", "/",
                                "100kb", "/",
                                "100kb", ".", numchr, ".annot.gz")))
  df2 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "Robospan_mean", "/",
                                "ABC", "/",
                                "ABC", ".", numchr, ".annot.gz")))
  df3 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "pRobospan_Mean_Blood", "/",
                                "100kb", "/",
                                "100kb", ".", numchr, ".annot.gz")))
  df4 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "pRobospan_mean", "/",
                                "ABC", "/",
                                "ABC", ".", numchr, ".annot.gz")))
  df5 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "pRobospan_mean", "/",
                                "Promoter_Gazal", "/",
                                "Promoter_Gazal", ".", numchr, ".annot.gz")))
  df6 = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES/",
                                "pRobospan_mean", "/",
                                "TSS", "/",
                                "TSS", ".", numchr, ".annot.gz")))

  newdf = cbind.data.frame(df1, df2[,5], df3[,5], df4[,5], df5[,5], df6[,5])
  colnames(newdf) = c(colnames(df1)[1:4], "Robospan_100kb", "Robospan_ABC",
                      "pRobospan_Blood_100kb", "pRobospan_ABC", "pRobospan_Promoter",
                      "pRobospan_TSS")

  if(!dir.exists(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/", "FS0"))){
    dir.create(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/", "FS0"))
  }
  write.table(newdf,
              file = gzfile(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/",
                                   "FS0", "/",
                                   "FS0", ".",
                                   numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  cat("We are at chrom:", numchr, "\n")

}

for(numchr in 1:22){
  df = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/",
                                "FS1", "/",
                                "FS1", ".", numchr, ".annot.gz")))
  newdf = df[,-c(5)]
  if(!dir.exists(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/", "FS2"))){
    dir.create(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/", "FS2"))
  }
  write.table(newdf,
              file = gzfile(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/FigJoints/",
                                   "FS2", "/",
                                   "FS2", ".",
                                   numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  cat("We are at chrom:", numchr, "\n")
}
