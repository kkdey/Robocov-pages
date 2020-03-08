
robocov_predict = get(load("/Users/kushaldey/Documents/Robocov-pages/output/Robocov_all_genes.rda"))
mean_vec = c()
for(gene in 1:length(dimnames(robocov_predict)[[3]])){
  brain_mat = robocov_predict[8:20, 8:20, gene]
  complement_mat = robocov_predict[-(8:20), -(8:20), gene]
  mean_vec = c(mean_vec, median(complement_mat[lower.tri(complement_mat)]) - median(brain_mat[lower.tri(brain_mat)]))
}
names(mean_vec) = dimnames(robocov_predict)[[3]]

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id",
                attributes= c("ensembl_gene_id","hgnc_symbol"),
                values=dimnames(robocov_predict)[[3]],mart= mart)
G_list = G_list[which(G_list[,2] != ""),]

idx = match(intersect(G_list[,1], dimnames(robocov_predict)[[3]]), dimnames(robocov_predict)[[3]])
mean_vec2 = mean_vec[idx]
names(mean_vec2) = G_list[match(intersect(G_list[,1], dimnames(robocov_predict)[[3]]), G_list[,1]), 2]

genes =  names(mean_vec2)[order(mean_vec2, decreasing=T)[1:2100]]
tmp = cbind.data.frame(genes, 1)
write.table(tmp, file = "/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/NonBrain_Robocov_mean_genes_top10.txt",
            col.names=F, row.names=F, sep = "\t", quote=F)

