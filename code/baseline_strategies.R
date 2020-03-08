ll = list.files("/n/groups/price/kushal/Mouse_Humans/data/ANNOTATIONS/GENE_SCORES2/SEG_GTEx")[-c(3,5,7)]
for(numchr in 1:22){
  base = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/Baselines/",
                                 "baselineLD_v2.1_S2G_Oct9", "/",
                                 "baselineLD", ".", numchr, ".annot.gz")))
  annot=c()
  temp=c()
  for(num in 1:length(ll)){
    df = data.frame(fread(paste0("zcat /n/groups/price/kushal/Mouse_Humans/data/ANNOTATIONS/GENE_SCORES2/",
                                 "SEG_GTEx",  "/",
                                 ll[num], "/", ll[num], ".", numchr, ".annot.gz")))
    annot = cbind(annot, df[,5])
  }
  colnames(annot) = paste0("SEG_", ll)

  base2 = cbind.data.frame(base, annot)

  if(!dir.exists(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/Baselines/", "baselineLD_v2.1_S2G_SEG_Nov23"))){
    dir.create(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/Baselines/", "baselineLD_v2.1_S2G_SEG_Nov23"))
  }
  write.table(base2,
              file = gzfile(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/Baselines/",
                                   "baselineLD_v2.1_S2G_SEG_Nov23", "/",
                                   "baselineLD", ".",
                                   numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)

  cat("We are at chrom:", numchr, "\n")
}


for(numchr in 1:22){
  base = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/Baselines/",
                                 "baselineLD_v2.1_S2G_Oct9", "/",
                                 "baselineLD", ".", numchr, ".annot.gz")))
  base2 = base[,c(1:92,96)]
  if(!dir.exists(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/Baselines/", "baselineLD_v2.1_S2G_Nov26"))){
    dir.create(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/Baselines/", "baselineLD_v2.1_S2G_Nov26"))
  }
  write.table(base2,
              file = gzfile(paste0("/n/groups/price/kushal/Robocov/data/ANNOTATIONS/Baselines/",
                                   "baselineLD_v2.1_S2G_Nov26", "/",
                                   "baselineLD", ".",
                                   numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)

  cat("We are at chrom:", numchr, "\n")
}


