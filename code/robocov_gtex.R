

robocov_predict = get(load("/Users/kushaldey/Documents/Robocov-pages/output/Cor_pairwise_all_genes.rda"))

library(corrplot)
library(ggplot2)
corrplot(robocov_predict[,,"ENSG00000244734"],  diag = TRUE,
         col = colorRampPalette(c("blue", "white", "red"))(200),
         tl.pos = "td", tl.cex = 0.4, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")

mean_vec = c()
num_vec = c()
for(mm in 1:length(dimnames(robocov_predict)[[3]])){
  tmp_mat = robocov_predict[,,mm]
  mean_vec = c(mean_vec, mean(robocov_predict[,,mm]))
  num_vec = c(num_vec,
              length(which(tmp_mat[lower.tri(tmp_mat)] > 0.1))/length(tmp_mat[lower.tri(tmp_mat)]))
  cat("We are at gene:", mm, "\n")
}

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id",
                attributes= c("ensembl_gene_id","hgnc_symbol"),
                values=dimnames(robocov_predict)[[3]],mart= mart)
G_list = G_list[which(G_list[,2] != ""),]

write.table(G_list, "/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/ensembl_and_hgnc_symbols.txt",
            row.names=F, col.names=F, quote=F, sep = "\t")

idx = match(intersect(G_list[,1], dimnames(robocov_predict)[[3]]), dimnames(robocov_predict)[[3]])
mean_vec2 = mean_vec[idx]
names(mean_vec2) = G_list[match(intersect(G_list[,1], dimnames(robocov_predict)[[3]]), G_list[,1]), 2]

num_vec2 = num_vec[idx]
names(num_vec2) = G_list[match(intersect(G_list[,1], dimnames(robocov_predict)[[3]]), G_list[,1]), 2]

genes2 = names(num_vec2)[order(num_vec2, decreasing=T)[1:2100]]
tmp = cbind.data.frame(genes2, 1)
write.table(tmp, file = "/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/Cor_num_genes_top10.txt",
            col.names=F, row.names=F, sep = "\t", quote=F)

genes =  names(mean_vec2)[order(mean_vec2, decreasing=T)[1:2100]]
tmp = cbind.data.frame(genes, 1)
write.table(tmp, file = "/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/Cor_mean_genes_top10.txt",
            col.names=F, row.names=F, sep = "\t", quote=F)




#tmp1 = read.table("/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/Cor_mean_genes_top10.txt")
#tmp2 = read.table("/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/Robocov_mean_genes_top10.txt")









