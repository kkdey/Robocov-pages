

cor_gtex = get(load("/Users/kushaldey/Documents/Robocov-pages/data/Cor_pairwise_all_genes.rda"))

corspan = apply(cor_gtex, 3, sum)/(53*53)
plot(density(corspan), xlab = "Corspan score", ylab = "density")

ensembl_gene_symbol = read.table("/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/ensembl_and_hgnc_symbols.txt")
head(ensembl_gene_symbol)
gene_symbols = ensembl_gene_symbol[match(names(corspan), ensembl_gene_symbol[,1]), 2]
corspan2 = corspan[which(!is.na(gene_symbols))]
names(corspan2) = gene_symbols[which(!is.na(gene_symbols))]

genes = names(corspan2)[order(corspan2, decreasing = T)[1:1600]]
df = cbind.data.frame(genes, 1)
write.table(df, file = "/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/Corspan_mean.txt",
            row.names = F, col.names = F, sep = "\t", quote=F)


corspan = apply(cor_gtex[,53,], 2, mean)
gene_symbols = ensembl_gene_symbol[match(names(corspan), ensembl_gene_symbol[,1]), 2]
corspan2 = corspan[which(!is.na(gene_symbols))]
names(corspan2) = gene_symbols[which(!is.na(gene_symbols))]
genes = names(corspan2)[order(corspan2, decreasing = T)[1:1600]]
df = cbind.data.frame(genes, 1)
write.table(df, file = "/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/Corspan_Mean_Blood.txt",
            row.names = F, col.names = F, sep = "\t", quote=F)

