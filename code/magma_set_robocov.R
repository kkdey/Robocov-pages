

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
robocov_gtex = get(load("/Users/kushaldey/Documents/Robocov-pages/data/Cor_pairwise_all_genes.rda"))

ensembl.genes = dimnames(robocov_gtex)[[3]]

genes_df <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=ensembl.genes,
  mart=mart)

unique_names = unique(genes_df[,1])
genes_df2 = genes_df[match(unique_names, genes_df[,1]), ]

ll= list.files("~/Documents/Robocov-pages/output/NMF_run_5_probocov/", pattern = ".txt")
outll = list()
for(numl in 1:length(ll)){
  genes =as.character(read.delim(paste0("~/Documents/Robocov-pages/output/NMF_run_5_probocov/", ll[numl]), header = F)[,1])
  entrez = genes_df2[match(genes, genes_df2[,1]), 2]
  entrez = entrez[!is.na(entrez)]
  outll[[numl]] = entrez
}
names(outll) = paste0("probocov_nmf_", ll)
outll1 = outll



ll= list.files("~/Documents/Robocov-pages/output/NMF_run_5_robocov/", pattern = ".txt")
outll = list()
for(numl in 1:length(ll)){
  genes =as.character(read.delim(paste0("~/Documents/Robocov-pages/output/NMF_run_5_robocov/", ll[numl]), header = F)[,1])
  entrez = genes_df2[match(genes, genes_df2[,1]), 2]
  entrez = entrez[!is.na(entrez)]
  outll[[numl]] = entrez
}
names(outll) = paste0("robocov_nmf_", ll)
outll2 = outll


ll= list.files("~/Documents/Robocov-pages/output/NMF_run_5_standard/", pattern = ".txt")
outll = list()
for(numl in 1:length(ll)){
  genes =as.character(read.delim(paste0("~/Documents/Robocov-pages/output/NMF_run_5_standard/", ll[numl]), header = F)[,1])
  entrez = genes_df2[match(genes, genes_df2[,1]), 2]
  entrez = entrez[!is.na(entrez)]
  outll[[numl]] = entrez
}
names(outll) = paste0("standard_nmf_", ll)
outll3 = outll

total_ll = c(outll1, outll2, outll3)

tt = c()
for(mm in 1:length(total_ll)){
  tt = c(tt, paste0(c(names(total_ll)[mm],
                      total_ll[[mm]][!is.na(total_ll[[mm]])]), collapse = "\t"))
}

writeLines(tt, "~/Documents/COVID19/data/MAGMA/magma_sets/robocov_nmf.sets",
           sep = "\n")





ll= list.files("~/Documents/Robocov-pages/output/Cor_hub_specific/", pattern = ".txt")
outll = list()
for(numl in 1:length(ll)){
  genes =as.character(read.delim(paste0("~/Documents/Robocov-pages/output/Cor_hub_specific/", ll[numl]), header = F)[,1])
  entrez = genes_df2[match(genes, genes_df2[,1]), 2]
  entrez = entrez[!is.na(entrez)]
  outll[[numl]] = entrez
}
names(outll) = paste0("cor_specific_", ll)
outll4 = outll

ll= list.files("~/Documents/Robocov-pages/output/Robocov_hub_specific/", pattern = ".txt")
outll = list()
for(numl in 1:length(ll)){
  genes =as.character(read.delim(paste0("~/Documents/Robocov-pages/output/Robocov_hub_specific/", ll[numl]), header = F)[,1])
  entrez = genes_df2[match(genes, genes_df2[,1]), 2]
  entrez = entrez[!is.na(entrez)]
  outll[[numl]] = entrez
}
names(outll) = paste0("robocov_specific_", ll)
outll5 = outll

ll= list.files("~/Documents/Robocov-pages/output/Robocov_precision_hub_specific/", pattern = ".txt")
outll = list()
for(numl in 1:length(ll)){
  genes =as.character(read.delim(paste0("~/Documents/Robocov-pages/output/Robocov_precision_hub_specific/", ll[numl]), header = F)[,1])
  entrez = genes_df2[match(genes, genes_df2[,1]), 2]
  entrez = entrez[!is.na(entrez)]
  outll[[numl]] = entrez
}
names(outll) = paste0("probocov_specific_", ll)
outll6 = outll


total_ll = c(outll4, outll5, outll6)

tt = c()
for(mm in 1:length(total_ll)){
  tt = c(tt, paste0(c(names(total_ll)[mm],
                      total_ll[[mm]][!is.na(total_ll[[mm]])]), collapse = "\t"))
}

writeLines(tt, "~/Documents/COVID19/data/MAGMA/magma_sets/robocov_hub_specific.sets",
           sep = "\n")







