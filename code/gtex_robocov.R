options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
group = as.numeric(toString(args[1]))

library(CVXR)
library(Robocov)

dat = get(load("/n/groups/price/kushal/Robocov/data/person_tissue_genes_voom.rda"))
num_samples_per_tissue = apply(dat[,,1], 2, function(x) return(length(which(!is.na(x)))))
gene_names_gtex = as.character(read.table("/n/groups/price/kushal/Robocov/data/gene_names_GTEX_V6.txt")[,1])
gene_names_gtex = as.character(sapply(gene_names_gtex, function(z) return(strsplit(z, "[.]")[[1]][1])))




chunk <- function(x, n) split(x, sort(rank(x) %% n))
block = chunk(1:dim(dat)[[3]], 20)
ids = block[[group]]

robomat = array(0, c(dim(dat)[2], dim(dat)[2], length(ids)))
for(n in 1:length(ids)){
  out = suppressMessages(suppressWarnings(Robocov_box(data_with_missing = dat[,,ids[n]], loss = "elasticnet")))
  robomat[,,n] = out
  cat("We are at gene", n, "\n")
}

saver(robomat, file = paste0("/n/groups/price/kushal/Robocov/data/", "Robocov_", group, ".rda"))

## ====================================================================================================================###


##   Post-processing the results from running the above script  (typically takes 40 mins across 20 batches each)    #####

library(abind)
tot_cov = c()
for(mm in 1:20){
  tmp = get(load(paste0("/n/groups/price/kushal/Robocov/data/", "Robocov_", mm, ".rda")))
  tot_cov = abind(tot_cov, tmp, along=3)
  cat("We are at chrunk:", mm, "\n")
}

dimnames(tot_cov)[[3]] = gene_names_gtex
dimnames(tot_cov)[[2]] = dimnames(dat)[[2]]
dimnames(tot_cov)[[1]] = dimnames(dat)[[2]]
save(tot_cov, file = "/n/groups/price/kushal/Robocov/data/Robocov_Box_all_genes.rda")




