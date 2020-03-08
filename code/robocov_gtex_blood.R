
robocov_predict = get(load("/Users/kushaldey/Documents/Robocov-pages/output/Robocov_all_genes.rda"))
blood_mat = robocov_predict[53, -53, ]
mean_vec = colMeans(blood_mat)
num_vec = apply(blood_mat, 2, function(x){
  tt = length(which(abs(x) > 0.1))
  return(tt)
})

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id",
                attributes= c("ensembl_gene_id","hgnc_symbol"),
                values=dimnames(robocov_predict)[[3]],mart= mart)
G_list = G_list[which(G_list[,2] != ""),]

idx = match(intersect(G_list[,1], dimnames(robocov_predict)[[3]]), dimnames(robocov_predict)[[3]])
mean_vec2 = mean_vec[idx]
names(mean_vec2) = G_list[match(intersect(G_list[,1], dimnames(robocov_predict)[[3]]), G_list[,1]), 2]

num_vec2 = num_vec[idx]
names(num_vec2) = G_list[match(intersect(G_list[,1], dimnames(robocov_predict)[[3]]), G_list[,1]), 2]


genes2 = names(num_vec2)[order(num_vec2, decreasing=T)[1:2100]]
tmp = cbind.data.frame(genes2, 1)
write.table(tmp, file = "/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/Blood_Robocov_num_genes_top10.txt",
            col.names=F, row.names=F, sep = "\t", quote=F)

genes =  names(mean_vec2)[order(mean_vec2, decreasing=T)[1:2100]]
tmp = cbind.data.frame(genes, 1)
write.table(tmp, file = "/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/Blood_Robocov_mean_genes_top10.txt",
            col.names=F, row.names=F, sep = "\t", quote=F)


