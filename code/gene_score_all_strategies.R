

###############################       Gene Score all strategies            ####################################

library(data.table)
annotation_cell = "/n/groups/price/kushal/Robocov/data/ANNOTATIONS/GENE_SCORES"
score_file = "/n/groups/price/kushal/Robocov/data/Gene_Scores/Hub_brain_robocov.txt"
annotation_folder = "Hub_brain_robocov"
bed_cell = "/n/groups/price/kushal/Robocov/data/BEDFILES"
bed_folder = "Hub_brain_robocov"


# gene_scores = read.delim(score_file, header=F)
#
# if(!dir.exists(paste0(annotation_cell, "/", annotation_folder))){
#   dir.create(paste0(annotation_cell, "/", annotation_folder))
# }
#
# for(numchr in 1:22){
#   cadd_scores = read.table(paste0("/n/groups/price/kushal/LDSC-Average/data/CADD_Scores/",
#                                   "CADD_ALL.", numchr, ".annot.gz"), header=T)
#   base = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/Baselines/",
#                                  "baselineLD_v2.1", "/", "baselineLD.", numchr, ".annot.gz")))
#   idx = match(cadd_scores$GeneName,gene_scores[,1])
#   binary_annot = rep(0, nrow(cadd_scores))
#   binary_annot[which(!is.na(idx))] = gene_scores[idx[!is.na(idx)],2]
#
#   temp = binary_annot[match(base$BP, cadd_scores$Pos)]
#   temp[is.na(temp)] = 0
#
#   if(!dir.exists(paste0(annotation_cell, "/", annotation_folder, "/", "GENE_5KB"))){
#     dir.create(paste0(annotation_cell, "/", annotation_folder, "/", "GENE_5KB"))
#   }
#
#   newdf = cbind.data.frame(base[,1:4], temp)
#   colnames(newdf) = c(colnames(base)[1:4], "GENE_5KB")
#
#   write.table(newdf,
#               file = gzfile(paste0(annotation_cell, "/", annotation_folder, "/",
#                                    "GENE_5KB/",
#                                    "GENE_5KB.", numchr, ".annot.gz")),
#               quote=FALSE, row.names=FALSE)
#
#   cat("We are at chr:", numchr, "\n")
# }
#
# for(numchr in 1:22){
#   df = data.frame(fread(paste0("zcat ", annotation_cell, "/", annotation_folder, "/",
#                                "GENE_5KB/", "GENE_5KB", ".", numchr, ".annot.gz")))
#   base = data.frame(fread(paste0("zcat /n/groups/price/kushal/Robocov/data/ANNOTATIONS/Baselines/",
#                                  "baselineLD_v2.1", "/", "baselineLD.", numchr, ".annot.gz")))
#   newdf1 = cbind.data.frame(df[,1:4], df[,5]*base$TSS_Hoffman)
#   colnames(newdf1) = c(colnames(base)[1:4], "AN")
#   newdf2 = cbind.data.frame(df[,1:4], df[,5]*base$Promoter_UCSC)
#   colnames(newdf2) = c(colnames(base)[1:4], "AN")
#   newdf3 = cbind.data.frame(df[,1:4], df[,5]*base$Coding_UCSC)
#   colnames(newdf3) = c(colnames(base)[1:4], "AN")
#
#   if(!dir.exists(paste0(annotation_cell, "/", annotation_folder, "/", "TSS"))){
#     dir.create(paste0(annotation_cell, "/", annotation_folder, "/", "TSS"))
#   }
#
#   write.table(newdf1,
#               file = gzfile(paste0(annotation_cell, "/", annotation_folder, "/",
#                                    "TSS/",
#                                    "TSS.", numchr, ".annot.gz")),
#               quote=FALSE, row.names=FALSE)
#
#   if(!dir.exists(paste0(annotation_cell, "/", annotation_folder, "/", "PROMOTER"))){
#     dir.create(paste0(annotation_cell, "/", annotation_folder, "/", "PROMOTER"))
#   }
#
#   write.table(newdf2,
#               file = gzfile(paste0(annotation_cell, "/", annotation_folder, "/",
#                                    "PROMOTER/",
#                                    "PROMOTER.", numchr, ".annot.gz")),
#               quote=FALSE, row.names=FALSE)
#
#   if(!dir.exists(paste0(annotation_cell, "/", annotation_folder, "/", "CODING"))){
#     dir.create(paste0(annotation_cell, "/", annotation_folder, "/", "CODING"))
#   }
#
#   write.table(newdf3,
#               file = gzfile(paste0(annotation_cell, "/", annotation_folder, "/",
#                                    "CODING/",
#                                    "CODING.", numchr, ".annot.gz")),
#               quote=FALSE, row.names=FALSE)
#
#   cat("We are at chr:", numchr, "\n")
#
# }


# source("/n/groups/price/kushal/Mouse_Humans/code/ALL_bedgraph_calc.R")
#
# out1 = ABC_bedgraph_calc(scores,
#                          output_cell = paste0(bed_cell, "/", bed_folder),
#                          output_bed = "ABC.bed")
#
# out2 = PCHiC_bedgraph_calc(scores,
#                            output_cell = paste0(bed_cell, "/", bed_folder),
#                            output_bed = "PCHiC.bed")
#
# out3 = Promoter_Gazal_bedgraph_calc(scores,
#                                     output_cell = paste0(bed_cell, "/", bed_folder),
#                                     output_bed = "Promoter_Gazal.bed")
#
# out4 = Roadmap_Enhancer_bedgraph_calc(scores,
#                                       output_cell = paste0(bed_cell, "/", bed_folder),
#                                       output_bed = "Roadmap_Enhancer.bed")
#
# out5 = Yoshida_bedgraph_calc(scores,
#                              output_cell = paste0(bed_cell, "/", bed_folder),
#                              output_bed = "Yoshida.bed")
#
# out6 = eQTLCPP_bedgraph_calc(scores,
#                              output_cell = paste0(bed_cell, "/", bed_folder),
#                              output_bed = "FinemapBloodeQTL.bed")


if(!dir.exists(paste0(bed_cell, "/", bed_folder))){
  dir.create(paste0(bed_cell, "/", bed_folder))
}

gene_scores = read.delim(score_file, header=F)

scores = gene_scores[,2]
names(scores) = gene_scores[,1]

df = read.table("/n/groups/price/kushal/Mouse_Humans/data/Gene_10kb.txt", header=T)
df[which(df[,2] < 0), 2] = 0
score = gene_scores[match(df$gene, gene_scores[,1]), 2]
score[is.na(score)] = 0

temp = cbind.data.frame(df[,1:3], score)
temp2 = temp[which(temp$score != 0),]
write.table(temp2, file = paste0(bed_cell, "/", bed_folder, "/", "10kb.bed"),
            quote=F, sep = "\t", row.names=F, col.names=F)


df = read.table("/n/groups/price/kushal/Mouse_Humans/data/Gene_100kb.txt", header=T)
df[which(df[,2] < 0), 2] = 0
score = gene_scores[match(df$gene, gene_scores[,1]), 2]
score[is.na(score)] = 0

temp = cbind.data.frame(df[,1:3], score)
temp2 = temp[which(temp$score != 0),]
write.table(temp2, file = paste0(bed_cell, "/", bed_folder, "/", "100kb.bed"),
            quote=F, sep = "\t", row.names=F, col.names=F)







