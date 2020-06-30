

robocov_gtex = get(load("/Users/kushaldey/Documents/Robocov-pages/data/Robocov_Precision_all_genes.rda"))
cell = "Robocov_precision_hub_specific"
num_genes=500

tissue_gtex = robocov_gtex[8:20, 8:20, ]
nontissue_gtex = robocov_gtex[-(8:20), -(8:20), ]
tt = apply(tissue_gtex, 3, mean) - apply(nontissue_gtex, 3, mean)
tt2 = tt[order(tt, decreasing = T)[1:num_genes]]
matt3 = apply(robocov_gtex[,,order(tt, decreasing = T)[1:num_genes]], c(1,2), mean)
corrplot(matt3,   diag = TRUE,
         col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
         tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "lower", tl.srt=45)
write.table(names(tt2), file = paste0("~/Documents/Robocov-pages/output/", cell, "/",
                                      "brain_specific.txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)



tissue_gtex = robocov_gtex[32:33, 32:33, ]
nontissue_gtex = robocov_gtex[-(32:33), -(32:33), ]
tt = apply(tissue_gtex, 3, mean) - apply(nontissue_gtex, 3, mean)
tt2 = tt[order(tt, decreasing = T)[1:num_genes]]
matt3 = apply(robocov_gtex[,,order(tt, decreasing = T)[1:num_genes]], c(1,2), mean)
corrplot(matt3,   diag = TRUE,
         col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
         tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "lower", tl.srt=45)
write.table(names(tt2), file = paste0("~/Documents/Robocov-pages/output/", cell, "/",
                                      "heart_specific.txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)



tissue_gtex = robocov_gtex[1:2, 1:2, ]
nontissue_gtex = robocov_gtex[-(1:2), -(1:2), ]
tt = apply(tissue_gtex, 3, mean) - apply(nontissue_gtex, 3, mean)
tt2 = tt[order(tt, decreasing = T)[1:num_genes]]
matt3 = apply(robocov_gtex[,,order(tt, decreasing = T)[1:num_genes]], c(1,2), mean)
corrplot(matt3,   diag = TRUE,
         col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
         tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "lower", tl.srt=45)
write.table(names(tt2), file = paste0("~/Documents/Robocov-pages/output/", cell, "/",
                                      "adipose_specific.txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)


tissue_gtex = robocov_gtex[c(28:30), c(28:30), ]
nontissue_gtex = robocov_gtex[-c(28:30), -c(28:30), ]
tt = apply(tissue_gtex, 3, mean) - apply(nontissue_gtex, 3, mean)
tt2 = tt[order(tt, decreasing = T)[1:num_genes]]
matt3 = apply(robocov_gtex[,,order(tt, decreasing = T)[1:num_genes]], c(1,2), mean)
corrplot(matt3,   diag = TRUE,
         col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
         tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "lower", tl.srt=45)
write.table(names(tt2), file = paste0("~/Documents/Robocov-pages/output/", cell, "/",
                                      "esophagus_specific.txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)


tissue_gtex = robocov_gtex[c(26, 27, 46, 48), c(26, 27, 46, 48), ]
nontissue_gtex = robocov_gtex[-c(26, 27, 46, 48), -c(26, 27, 46, 48), ]
tt = apply(tissue_gtex, 3, mean) - apply(nontissue_gtex, 3, mean)
tt2 = tt[order(tt, decreasing = T)[1:num_genes]]
matt3 = apply(robocov_gtex[,,order(tt, decreasing = T)[1:num_genes]], c(1,2), mean)
corrplot(matt3,   diag = TRUE,
         col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
         tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "lower", tl.srt=45)
write.table(names(tt2), file = paste0("~/Documents/Robocov-pages/output/", cell, "/",
                                      "digestive_specific.txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)


tissue_gtex = robocov_gtex[c(4:6), c(4:6), ]
nontissue_gtex = robocov_gtex[-c(4:6), -c(4:6), ]
tt = apply(tissue_gtex, 3, mean) - apply(nontissue_gtex, 3, mean)
tt2 = tt[order(tt, decreasing = T)[1:num_genes]]
matt3 = apply(robocov_gtex[,,order(tt, decreasing = T)[1:num_genes]], c(1,2), mean)
corrplot(matt3,   diag = TRUE,
         col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
         tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "lower", tl.srt=45)
write.table(names(tt2), file = paste0("~/Documents/Robocov-pages/output/", cell, "/",
                                      "artery_specific.txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)


tissue_gtex = robocov_gtex[c(44:45), c(44:45), ]
nontissue_gtex = robocov_gtex[-c(44:45), -c(44:45), ]
tt = apply(tissue_gtex, 3, mean) - apply(nontissue_gtex, 3, mean)
tt2 = tt[order(tt, decreasing = T)[1:num_genes]]
matt3 = apply(robocov_gtex[,,order(tt, decreasing = T)[1:num_genes]], c(1,2), mean)
corrplot(matt3,   diag = TRUE,
         col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
         tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "lower", tl.srt=45)
write.table(names(tt2), file = paste0("~/Documents/Robocov-pages/output/", cell, "/",
                                      "artery_specific.txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)


tissue_gtex = robocov_gtex[-(8:20), -(8:20), ]
nontissue_gtex = robocov_gtex[c(8:20), c(8:20), ]
tt = apply(tissue_gtex, 3, mean) - apply(nontissue_gtex, 3, mean)
tt2 = tt[order(tt, decreasing = T)[1:num_genes]]
matt3 = apply(robocov_gtex[,,order(tt, decreasing = T)[1:num_genes]], c(1,2), mean)
corrplot(matt3,   diag = TRUE,
         col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
         tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "lower", tl.srt=45)
write.table(names(tt2), file = paste0("~/Documents/Robocov-pages/output/", cell, "/",
                                      "nonbrain_specific.txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)

tissue_gtex = robocov_gtex[c(51, 52), c(51, 52), ]
nontissue_gtex = robocov_gtex[-c(51, 52), -c(51, 52), ]
tt = apply(tissue_gtex, 3, mean) - apply(nontissue_gtex, 3, mean)
tt2 = tt[order(tt, decreasing = T)[1:num_genes]]
matt3 = apply(robocov_gtex[,,order(tt, decreasing = T)[1:num_genes]], c(1,2), mean)
corrplot(matt3,   diag = TRUE,
         col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
         tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "lower", tl.srt=45)
write.table(names(tt2), file = paste0("~/Documents/Robocov-pages/output/", cell, "/",
                                      "sex_specific_female.txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)

tissue_gtex = robocov_gtex[c(7, 43, 49), c(7, 43, 49), ]
nontissue_gtex = robocov_gtex[-c(7, 43, 49), -c(7, 43, 49), ]
tt = apply(tissue_gtex, 3, mean) - apply(nontissue_gtex, 3, mean)
tt2 = tt[order(tt, decreasing = T)[1:num_genes]]
matt3 = apply(robocov_gtex[,,order(tt, decreasing = T)[1:num_genes]], c(1,2), mean)
corrplot(matt3,   diag = TRUE,
         col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
         tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "lower", tl.srt=45)
write.table(names(tt2), file = paste0("~/Documents/Robocov-pages/output/", cell, "/",
                                      "sex_specific_male.txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)

tissue_gtex = robocov_gtex
tt = apply(tissue_gtex, 3, mean)
tt2 = tt[order(tt, decreasing = T)[1:num_genes]]
matt3 = apply(robocov_gtex[,,order(tt, decreasing = T)[1:num_genes]], c(1,2), mean)
corrplot(matt3,   diag = TRUE,
         col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
         tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "lower", tl.srt=45)
write.table(names(tt2), file = paste0("~/Documents/Robocov-pages/output/", cell, "/",
                                      "all_specific.txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)


tissue_gtex = robocov_gtex
tt = apply(tissue_gtex, 3, mean)
tt2 = tt[order(tt, decreasing = F)[1:num_genes]]
matt3 = apply(robocov_gtex[,,order(tt, decreasing =F)[1:num_genes]], c(1,2), mean)
corrplot(matt3,   diag = TRUE,
         col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
         tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "lower", tl.srt=45)
write.table(names(tt2), file = paste0("~/Documents/Robocov-pages/output/", cell, "/",
                                      "none_specific.txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)






corrplot(robocov_gtex[,,"ENSG00000167768"],   diag = TRUE,
         col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
         tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "lower", tl.srt=45)
