
library(NNLM)
robocov_gtex = get(load("/Users/kushaldey/Documents/Robocov-pages/data/Cor_pairwise_all_genes.rda"))

robocov_mat = apply(robocov_gtex, 3, function(x){
  y = x[lower.tri(x)]
  return(y)
})

robocov_mat[robocov_mat < 0] = 0
nn = NNLM::nnmf(robocov_mat, k = 5)

cor_fac = cor(t(nn$H))
vec = c()
for(mm in 2:ncol(cor_fac)){
  xx = cor_fac[1:mm, 1:mm]
  if(max(xx[lower.tri(xx)]) < 0.5){
    vec = c(vec, mm)
  }
}


cor_fac[cor_fac > 0.5] = 0
corrplot(cor_fac,   diag = TRUE,
         col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
         tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "lower", tl.srt=45)

factors = nn$H
factors = t(apply(factors, 1, function(x) return(x/max(x))))

for(qq in 1:nrow(factors)){
  gene_indices = which(factors[qq,] > 0.5)
  matt2 = robocov_gtex[,,gene_indices]
  write.table(names(gene_indices),
              file = paste0("~/Documents/Robocov-pages/output/NMF_run_5_standard/NMF_run_5_clus_",
              qq, ".txt"), row.names = F, col.names = F, sep = "\t", quote=F)
  matt3 = apply(matt2, c(1,2), mean)
  corrplot(matt3,   diag = TRUE,
           col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
           tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
           rect.col = "white",na.label.col = "white",
           method = "color", type = "lower", tl.srt=45)
  cat("We are at factor:", qq, "\n")
}

target_loss = c()
for(kk in 1:20){
  nn = NNLM::nnmf(robocov_mat, k = kk)
  target_loss = c(target_loss, tail(nn$target.loss, 1))
}






