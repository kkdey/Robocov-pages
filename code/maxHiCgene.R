

########################  MaxHiCGene scores   ###########################

gene_summaries = read.table(gzfile(paste0("/n/groups/price/kushal/HiCDeep/data/", "gene_summary_hiC.gz")),
                            header = T)
files_to_include = rownames(gene_summaries)
tss_summaries = gene_summaries$tss

annot1 = c()
annot2 = c()
annot3 = c()

for(numfile in 1:nrow(gene_summaries)){
  hiC_scores_data = data.frame(data.table::fread(paste0("zcat /n/groups/price/steven/h2gene/DATA/HiC_Engreitz/",
                                                        files_to_include[numfile])))
  hiC_scores_data_2 = hiC_scores_data
  indices1 = which(hiC_scores_data[,2] > (tss_summaries[numfile] - 1e04) & hiC_scores_data[,2] < (tss_summaries[numfile] + 1e04))
  hiC_scores_data_2[indices1, 4] = 0
  annot1 = c(annot1, mean(sort(hiC_scores_data_2[,4], decreasing=TRUE)[1:3]))

  hiC_scores_data_3 = hiC_scores_data
  indices2 = which(hiC_scores_data[,2] > (tss_summaries[numfile] - 1e05) & hiC_scores_data[,2] < (tss_summaries[numfile] + 1e05))
  hiC_scores_data_3[indices2, 4] = 0
  annot2 = c(annot2, mean(sort(hiC_scores_data_3[,4], decreasing=TRUE)[1:3]))

  hiC_scores_data_4 = hiC_scores_data
  indices3 = which(hiC_scores_data[,2] > (tss_summaries[numfile] - 1e06) & hiC_scores_data[,2] < (tss_summaries[numfile] + 1e06))
  hiC_scores_data_4[indices3, 4] = 0
  annot3 = c(annot3, mean(sort(hiC_scores_data_4[,4], decreasing=TRUE)[1:3]))

  cat("We recorded data of percentage :", round(100*(numfile/length(files_to_include))), "\n")
}

df = cbind.data.frame(gene_summaries, annot1, annot2, annot3)
colnames(df) = c(colnames(gene_summaries), "TSS-1E04", "TSS-1E05", "TSS-1E06")

write.table(df,
            file = gzfile(paste0("/n/groups/price/kushal/HiCDeep/data/",
                                 "MaxHiC_GENE_Scores.gz")),
            quote=FALSE, row.names=TRUE)


